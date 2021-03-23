"""Calculate the ranking of the hits."""

import os
import json
import requests
import logging
import pandas as pd

from datetime import datetime
from tqdm import tqdm
from json.decoder import JSONDecodeError
from ebel_rest import query as rest_query

from dif.constants import PUBCHEM_BIOASSAY_API, UNIPROT_ID, CT_MAPPER, ACTION_MAPPER, CLINICAL_TRIAL_FROM_DRUG, \
    ASSOCIATED_PATHWAYS
from dif.defaults import BIOASSAY_CACHE
from dif.finder import InteractorFinder

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Ranker:
    """Class for taking an InteractorFinder object with saved results and annotating said results with ranking
    metadata."""

    def __init__(self, finder: InteractorFinder, penalty: int = -1, reward: int = 1):
        """Should be initialized with an InteractorFinder object that has results saved."""
        self.finder = finder
        self.penalty = penalty
        self.reward = reward

        if self.finder.results is None:
            raise ValueError("InteractorFinder.results is empty! First use InteractorFinder().find_interactors() or "
                             "InteractorFinder().druggable_interactors().")

        self.results = finder.results
        self.interactors = finder.unique_interactors()
        self.drugs = finder.unique_drugs()

    @staticmethod
    def __check_bioassay_cache() -> dict:
        """Checks if BioAssay cache file present. If not, creates it. Returns current content."""
        if not os.path.exists(BIOASSAY_CACHE):
            open(BIOASSAY_CACHE, 'w').close()
            counts = dict()

        else:
            try:
                with open(BIOASSAY_CACHE, 'r') as bioassay_file:
                    counts = json.load(bioassay_file)

            except JSONDecodeError:  # File empty
                counts = dict()
                logger.info("Created new BioAssay cache file.")

        return counts

    def get_bioassay_counts(self) -> dict:
        """Queries the PubChem API to obtain the number of available BioAssays for each gene symbol in the list.

        Returns
        -------
        dict
            Key is the gene symbol, value is the number of BioAssays available.
        """
        counts = self.__check_bioassay_cache()
        symbols_missing_data = set(self.interactors) - set(counts.keys())

        if symbols_missing_data:
            for symbol in tqdm(symbols_missing_data, desc="Counting BioAssays for targets"):
                if symbol in counts:  # Already in counts
                    continue

                up_acc_results = rest_query.sql(UNIPROT_ID.format(symbol)).data

                if up_acc_results:
                    up_acc = up_acc_results[0]['uniprot.id']

                    filled = PUBCHEM_BIOASSAY_API.format(up_acc)
                    resp = requests.get(filled)
                    number_assays = len(resp.text.split("\n")) - 1  # Remove header
                    counts[symbol] = number_assays

        with open(BIOASSAY_CACHE, 'w') as bioassay_file:
            json.dump(counts, bioassay_file)
            logger.info(f"BioAssay data written to {BIOASSAY_CACHE}")

        return counts

    def score_patents(self, drugs_and_patents_raw: dict) -> dict:
        drug_patents_parsed = dict()
        for drug_name, patent_info in drugs_and_patents_raw.items():
            if patent_info is None:
                drug_patents_parsed[drug_name] = {'has_patent': False, 'expired': False,
                                                  'number': None, 'points': self.reward}

            else:
                patent_numbers = []
                expired_list = []
                patent_list = patent_info['patent'] if isinstance(patent_info['patent'], list) else [patent_info['patent']]
                for patent in patent_list:
                    patent_numbers.append(patent['number'])
                    expired_list.append(datetime.today() > datetime.strptime(patent['expires'], "%Y-%m-%d"))

                all_expired = all(expired_list)
                drug_patents_parsed[drug_name] = {'has_patent': bool(patent_numbers),
                                                  'expired': all_expired,
                                                  'patent_numbers': patent_numbers,
                                                  'points': self.reward if all_expired else self.penalty}

        return drug_patents_parsed

    def score_generics(self, drugs_and_product_info: dict) -> dict:
        """Gives scores to drugs in the hit list based on whether they have an approved generic version.

        Parameters
        ----------
        drugs_and_product_info: dict
            Keys are drug names and values are drug products.

        Returns
        -------
        generic_scores: dict
            Keys are drug names, values are generic results and point totals.
        """
        generic_scores = dict()
        for drug_name, product_info in drugs_and_product_info.items():
            has_generic = False
            has_approved_generic = False
            generic_info = dict()

            if product_info is None:
                pts = self.penalty

            else:
                if not isinstance(product_info, list):
                    product_info = [product_info]

                for product in product_info:
                    if product['generic'] == 'true':
                        has_generic = True
                        is_approved = True if product['approved'] == 'true' else False
                        product_data = {'is_approved': is_approved}
                        if is_approved:
                            has_approved_generic = True
                        generic_info[product['name']] = product_data

                if has_generic and has_approved_generic:
                    pts = self.reward
                else:
                    pts = self.penalty

            generic_scores[drug_name] = {'has_generic': has_generic, 'has_approved_generic': has_approved_generic,
                                         'generic_products': generic_info, 'points': pts}

        return generic_scores

    def score_homologs(self, tc_json: str, db_id_mapper: dict, tc_threshold: int = 0.95):
        """Produces a dictionary detailing structural homologs of every drugbank drug."""
        with open(tc_json, 'r') as mcsf:
            content = json.load(mcsf)
        homolog_scores = dict()
        for db_id, props in content.items():
            if 'similarities' in props and db_id in db_id_mapper:  # Check if in mapper else not in KG list
                drug_name = db_id_mapper[db_id]
                homolog_scores[drug_name] = {'homologs': {}, 'points': self.penalty}
                for db_id2, tc in props['similarities'].items():
                    if tc >= tc_threshold:
                        homolog_scores[drug_name]['homologs'][db_id2] = tc
                        homolog_scores[drug_name]['points'] = self.reward
        return homolog_scores

    # Drug/Target Relation and Interactor/pTau Relation Syngery
    @staticmethod
    def __check_contradictions(drug_rel_dict: dict) -> dict:
        """Check whether there are contradicting relation types between the interactor and pTAU."""
        tests = (('increases', 'decreases'),
                 ('increases', 'directly_decreases'),
                 ('decreases', 'directly_increases'),
                 ('directly_increases', 'directly_decreases'))

        for drug_pair, data in drug_rel_dict.items():
            rels = data['relation_type']
            contradictions = []
            for test in tests:
                contradiction_check = all(rel in rels for rel in test)
                contradictions.append(contradiction_check)

            drug_rel_dict[drug_pair]['contradictions'] = any(contradictions)

        return drug_rel_dict

    @staticmethod
    def __check_drug_action_contradictions(drug_actions: set) -> bool:
        """Checks if both 'positive_regulator' and 'negative_regulator' in the set of mapped drug/target relations."""
        contradiction = False
        if 'positive_regulator' and 'negative_regulator' in drug_actions:
            contradiction = True
        return contradiction

    @staticmethod
    def __compare_rels(drug_actions: set, int_rel: set) -> bool:
        """Returns True is drug/target relation and interactor/pTAU relation results in decrease of pTAU else False."""
        pos_drug, neg_drug = "positive_regulator", "negative_regulator"
        pos_int_rel = ['increases', 'directly_increases']
        neg_int_rel = ['decreases', 'directly_decreases']

        # drug -inhibits-> interactor -decreases-> pTAU: BAD
        if neg_drug in drug_actions and any(rel in int_rel for rel in neg_int_rel):
            return False

        # drug -inhibits-> interactor -increases-> pTAU: GOOD
        elif neg_drug in drug_actions and any(rel in int_rel for rel in pos_int_rel):
            return True

        # drug -activates-> interactor -increases-> pTAU: BAD
        elif pos_drug in drug_actions and any(rel in int_rel for rel in pos_int_rel):
            return False

        # drug -activates-> interactor -decreases-> pTAU: GOOD
        elif pos_drug in drug_actions and any(rel in int_rel for rel in neg_int_rel):
            return True

        else:
            print(drug_actions)
            print(int_rel)

    def score_drug_relationships(self, checked_drug_rels: dict) -> dict:
        """Adds the score to the drug/target pair based on whether it has information on the action (inhibitor/activator)
        and whether it makes sense based on the relationship of the interactor with pTAU.

        Parameters
        ----------
        checked_drug_rels: dict
            Keys are drug/target pairs, values are dictionaries realtion types.

        Returns
        -------
        checked_drug_rels: dict
            Keys are drug names, values are point values.
        """

        for drug_pair, data in checked_drug_rels.items():
            pts = 0
            if data['contradictions'] is True:
                pts += self.penalty
                checked_drug_rels[drug_pair]['points'] = pts  # Update drug_rels dict with points
                continue  # No point in continuing since drug/target action and interactor/pTau rel can't be compared

            rels = data['relation_type']  # imported as a set

            if data['actions'] is None:
                pts += self.penalty

            else:
                mapped_actions = {ACTION_MAPPER[rel] for rel in data['actions']}
                # First check if drug actions have contradictions
                if self.__check_drug_action_contradictions(mapped_actions) is False:
                    pts += self.penalty

                else:  # Compare synergy of drug actions with interactor/pTAU relations
                    if self.__compare_rels(mapped_actions, rels) is True:  # Good comparison
                        pts += self.reward

                    else:
                        pts += self.penalty

            checked_drug_rels[drug_pair]['points'] = pts  # Update drug_rels dict with points

        return checked_drug_rels

    @staticmethod
    def __collect_ct_info(db_id_mapper: dict) -> dict:
        """Collects clinical trial information for every identified drug.

        Parameters
        ----------
        db_id_mapper: dict
            Keys are DrugBank IDs and values are the drug names.

        Returns
        -------
        drug_ct_data: dict

        """
        drug_ct_data = dict()

        for db_id, drug_name in tqdm(db_id_mapper.items(), total=len(db_id_mapper), desc="Gathering clinical trial info"):
            results = rest_query.sql(CLINICAL_TRIAL_FROM_DRUG.format(db_id)).to_dict('records')
            if results:
                drug_ct_data[db_id] = {'drug_name': drug_name, 'ct_data': dict()}
                for r in results:
                    status = r['overall_status']
                    primary_condition = r['condition']  # List
                    mesh = r['mesh_conditions'] if 'mesh_conditions' in r else None  # List or None
                    conditions = set(primary_condition + mesh) if mesh is not None else set(primary_condition)
                    trial_drugs = r['drugs_in_trial'] if 'drugs_in_trial' in r else None
                    drug_ct_data[db_id]['ct_data'][r['trial_id']] = {'trial_status': status,
                                                                     'conditions': conditions,
                                                                     'drugs_in_trial': trial_drugs}
        return drug_ct_data

    def score_cts(self, db_id_mapper: dict, disease_keyword: str = "Alzheimer Disease",
                  similar_diseases: list = []) -> dict:
        """Scores drugs based on involvement in a clinical trial. Returned dictionary contains
        points and relevant AD-associated CT information."""

        ct_scores = dict()
        ct_info = self.__collect_ct_info(db_id_mapper)
        for db_id, metadata in ct_info.items():
            drug_name = metadata['drug_name']
            ct_scores[drug_name] = {'ct_data': dict()}
            pts = self.reward  # Default to reward unless changed

            if metadata['ct_data']:  # There is clinical trial data

                for trial_id, trial_data in metadata['ct_data'].items():
                    if disease_keyword in trial_data['conditions']:  # Disease-associated CT
                        ct_scores[drug_name]['ct_data'][trial_id] = trial_data  # Add data even if CT done
                        if trial_data['trial_status'] in CT_MAPPER['ongoing']:
                            pts = self.penalty

                    else:  # Not associated with primary disease
                        if set(trial_data['conditions']) & set(similar_diseases):  # Trial for similar disease
                            pts = self.reward * 2  # Double points
                        else:
                            if trial_data['trial_status'] in CT_MAPPER['ongoing']:
                                pts = self.penalty

            ct_scores[drug_name]['points'] = pts

        return ct_scores

    def count_associated_pathways(self) -> dict:
        """Goes through each interactor and checks if it is associated with a KEGG or Pathway Commons pathway."""
        pathways = dict()
        for protein in tqdm(self.interactors, desc="Checking pathways"):
            r = rest_query.sql(ASSOCIATED_PATHWAYS.format(protein, protein))
            count = r.iloc[0]['count']
            pathways[protein] = count

        return pathways

    def get_druggable_interactor_links(self) -> pd.DataFrame:
        """Collects information on all of the interactions for each druggable interactor."""
        pass
