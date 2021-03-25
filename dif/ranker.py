"""Calculate the ranking of the hits."""

import os
import json
import requests
import logging

from tqdm import tqdm
from typing import Optional
from datetime import datetime
from json.decoder import JSONDecodeError
from ebel_rest import query as rest_query

from dif.constants import *
from dif.utils import chunks
from dif.finder import InteractorFinder
from dif.defaults import BIOASSAY_CACHE, session
from dif.models import Trials, Edges, Patents, Products

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Ranker:
    """Class for taking an InteractorFinder object with saved results and annotating said results with ranking
    metadata."""

    def __init__(self, finder: InteractorFinder, penalty: int = -1, reward: int = 1):
        """Should be initialized with an InteractorFinder object that has results saved."""
        self.finder = finder
        self.__penalty = penalty
        self.__reward = reward

        if self.finder.results is None:
            raise ValueError("InteractorFinder.results is empty! First use InteractorFinder().find_interactors() or "
                             "InteractorFinder().druggable_interactors().")

        self.table = finder.results

        self.interactors = finder.unique_interactors()
        self.interactor_metadata = {int_name: dict() for int_name in self.interactors}

        self.drugs = finder.unique_drugs()
        self.drug_metadata = self.__compile_drug_metadata()
        self.drug_scores = self.__generate_ranking_dict()

    def rank(self):
        self.score_drugs()
        self.score_interactors()

    def score_drugs(self):
        """Wrapper method to parse raw drug metadata and calculate points for each ranking criteria."""
        self.score_drug_relationships()
        self.score_patents()
        self.score_generics()
        self.score_cts()
        # self.score_homologs()

    def score_interactors(self):
        """Wrapper method to parse raw interactor metadata and calculate points for each ranking criteria."""
        self.count_bioassays()
        self.count_edges()

    def __compile_drug_metadata(self) -> dict:
        metadata = self.__generate_ranking_dict()

        missing_pp_db_ids = set()
        for _, r in tqdm(self.table.iterrows(), total=len(self.table)):
            drug_name = r['drug']
            db_id = r['drugbank_id']
            interactor_name = r['interactor_name']
            metadata[drug_name][IDENTIFIERS] = {'drugbank_id': db_id,
                                                'chembl_id': r['chembl_id'],
                                                'pubchem_id': r['pubchem_id']}

            if r['relation_type'] == 'regulates':
                continue
            elif interactor_name in metadata[drug_name][INTERACTORS]:
                metadata[drug_name][INTERACTORS][interactor_name]['relation_type'].add(r['relation_type'])
            else:
                rel_data = {'relation_type': {r['relation_type']},  # Collect relation types in set for later comparison
                            'actions': r['drug_rel_actions'].split("|") if r['drug_rel_actions'] else None}
                metadata[drug_name][INTERACTORS][interactor_name] = rel_data

            patent_data = self.__query_db_patents(drug_name)
            product_data = self.__query_db_products(drug_name)

            if patent_data is None or product_data is None:
                missing_pp_db_ids.add(db_id)

            if patent_data is not None:
                metadata[drug_name][PATENTS] = patent_data

            if product_data is not None:
                metadata[drug_name][PRODUCTS] = product_data

        pp_data = self.__query_graphstore_patents_products(list(missing_pp_db_ids))
        for drug_name, pp_metadata in tqdm(pp_data.items()):
            metadata[drug_name][PATENTS] = pp_metadata['drug_patents']
            metadata[drug_name][PRODUCTS] = pp_metadata['drug_products']

        return metadata

    def __generate_ranking_dict(self) -> dict:
        """Compiles a generic empty dict to fill with raw and scored data."""
        empty_dict = {drug_name: {PATENTS: dict(),
                                  PRODUCTS: dict(),
                                  IDENTIFIERS: dict(),
                                  INTERACTORS: dict(),
                                  CLINICAL_TRIALS: dict()
                                  }
                      for drug_name in self.drugs}
        return empty_dict

    @staticmethod
    def __query_db_patents(drug_name: str) -> Optional[dict]:
        """Obtains patent information from SQLite DB for given gene drug name."""
        results = session().query(Patents).filter_by(drug_name=drug_name).all()
        assert len(results) < 2
        if results:
            hit = results[0]
            data = {'drug_name': hit.drug_name, 'has_patent': hit.has_patent,
                    'expired': hit.expired, 'patent_numbers': hit.patent_numbers.split("|")}
            return data

    @staticmethod
    def __query_db_products(drug_name: str) -> Optional[dict]:
        """Obtains product information from SQLite DB for given gene drug name."""
        results = session().query(Products).filter_by(drug_name=drug_name).all()
        assert len(results) < 2
        if results:
            hit = results[0]
            data = {'drug_name': hit.drug_name,
                    'has_generic': hit.has_generic,
                    'has_approved_generic': hit.has_approved_generic,
                    'generic_products': hit.generic_products.split("|")}
            return data

    @staticmethod
    def __query_graphstore_patents_products(db_ids: list) -> Optional[dict]:
        """Used to retrieve patent/product information for DrugBank IDs not in DB."""
        logger.info("Querying graphstore for patent/product information")
        sess = session()
        pp_data = dict()
        for id_list_chunk in chunks(db_ids, n=400):
            results = rest_query.sql(PATENTS_PRODUCTS.format(id_list_chunk)).data
            for hit in results:
                drug_name = hit.pop('name')
                pp_data[drug_name] = hit  # Remaining keys are drug_patents + drug_products

                # Add to DB
                patent_entry = {**hit['drug_patents'], **{'drug_name': drug_name}}
                product_entry = {**hit['drug_products'], **{'drug_name': drug_name}}

                sess.add(Patents(**patent_entry))
                sess.add(Products(**product_entry))

            sess.commit()

        return pp_data if pp_data else None

    def score_patents(self):
        drugs_and_patents_raw = {drug_name: dd[PATENTS] for drug_name, dd in self.drug_metadata.items()}
        for drug_name, patent_info in drugs_and_patents_raw.items():
            patent_numbers = []
            expired_list = []
            if patent_info:
                patent_list = patent_info['patent']

                if not isinstance(patent_info['patent'], list):
                    patent_list = [patent_info['patent']]

                for patent in patent_list:
                    patent_numbers.append(patent['number'])
                    expired_list.append(datetime.today() > datetime.strptime(patent['expires'], "%Y-%m-%d"))

                all_expired = all(expired_list)
                self.drug_scores[drug_name][PATENTS] = {'has_patent': bool(patent_numbers),
                                                        'expired': all_expired,
                                                        'patent_numbers': "|".join(patent_numbers),
                                                        POINTS: self.__reward if all_expired else self.__penalty}

            else:
                self.drug_scores[drug_name][PATENTS] = None

    def score_generics(self):
        """Gives scores to drugs in the hit list based on whether they have an approved generic version.

        Returns
        -------
        generic_scores: dict
            Keys are drug names, values are generic results and point totals.
        """
        drugs_and_product_info = {drug_name: dd[PRODUCTS] for drug_name, dd in self.drug_metadata.items()}
        for drug_name, product_info in drugs_and_product_info.items():
            has_generic = False
            has_approved_generic = False
            generic_info = dict()

            if product_info is None:
                pts = self.__penalty

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
                    pts = self.__reward
                else:
                    pts = self.__penalty

            self.drug_scores[drug_name][PRODUCTS] = {'has_generic': has_generic,
                                                     'has_approved_generic': has_approved_generic,
                                                     'generic_products': generic_info,
                                                     POINTS: pts}

    def score_homologs(self, tc_json: str, db_id_mapper: dict, tc_threshold: int = 0.95):
        """Produces a dictionary detailing structural homologs of every drugbank drug."""
        with open(tc_json, 'r') as mcsf:
            content = json.load(mcsf)
        homolog_scores = dict()
        for db_id, props in content.items():
            if 'similarities' in props and db_id in db_id_mapper:  # Check if in mapper else not in KG list
                drug_name = db_id_mapper[db_id]
                homolog_scores[drug_name] = {'homologs': {}, POINTS: self.__penalty}
                for db_id2, tc in props['similarities'].items():
                    if tc >= tc_threshold:
                        homolog_scores[drug_name]['homologs'][db_id2] = tc
                        homolog_scores[drug_name][POINTS] = self.__reward
        return homolog_scores

    # Drug/Target Relation and Interactor/pTau Relation Syngery
    @staticmethod
    def __check_target_interactor_contradictions(interactor_metadata: dict) -> bool:
        """Check whether there are contradicting relation types between the interactor and target."""
        tests = (('increases', 'decreases'),
                 ('increases', 'directly_decreases'),
                 ('decreases', 'directly_increases'),
                 ('directly_increases', 'directly_decreases'))

        relations = interactor_metadata['relation_type']
        contradictions = []
        for test in tests:
            contradiction_check = all(rel in relations for rel in test)
            contradictions.append(contradiction_check)

        return any(contradictions)

    @staticmethod
    def __check_drug_action_contradictions(drug_actions: set) -> bool:
        """Checks if both 'positive_regulator' and 'negative_regulator' in the set of mapped drug/target relations."""
        contradiction = False
        if 'positive_regulator' and 'negative_regulator' in drug_actions:
            contradiction = True
        return contradiction

    @staticmethod
    def __check_rel_synergy(drug_actions: set, int_rel: set) -> bool:
        """Returns True is drug/target relation and interactor/pTAU relation results in decrease of pTAU else False."""
        pos_drug, neg_drug = "positive_regulator", "negative_regulator"
        pos_int_rel = ['increases', 'directly_increases']
        neg_int_rel = ['decreases', 'directly_decreases']

        # TODO make generic so it works for any target
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

        else:  # drug_actions = 'neutral'
            return False

    def score_drug_relationships(self):
        """Adds the score to the drug/target pair based on whether it has information on the action
        (inhibitor/activator) and whether it makes sense based on the relationship of the interactor with target.
        Checks:
        * Whether there are contradicting edges between target/interactor (TIC) --> penalty
        * Whether there are contradicting drug actions between drug/interactor (DAC) --> penalty
        * If the drug action and target/interactor relationship synergize to produce desired effect (synergy)

        Returns
        -------
        checked_drug_rels: dict
            Keys are drug names, values are point values.
        """
        for drug_name, metadata in self.drug_metadata.items():
            for target_name, ti_metadata in metadata[INTERACTORS].items():
                self.drug_scores[drug_name][INTERACTORS][target_name] = {
                    TIC: False,
                    DAC: False,
                    'synergy': False,
                    POINTS: 0
                }
                tic = self.__check_target_interactor_contradictions(ti_metadata)  # Target/interactor

                pts = 0
                if tic is True:  # There are contradicting edges between target and interactor
                    pts += self.__penalty
                    self.drug_scores[drug_name][INTERACTORS][target_name][TIC] = tic
                    self.drug_scores[drug_name][INTERACTORS][target_name][POINTS] = pts
                    continue  # No point in continuing since drug/target action and int/target rel can't be compared

                rels = ti_metadata['relation_type']  # imported as a set

                if ti_metadata['actions'] is None:
                    pts += self.__penalty

                else:
                    mapped_actions = {ACTION_MAPPER[rel] for rel in ti_metadata['actions']}
                    # First check if drug actions have contradictions
                    dac = self.__check_drug_action_contradictions(mapped_actions)
                    if dac is True:
                        self.drug_scores[drug_name][INTERACTORS][target_name][DAC] = dac
                        pts += self.__penalty

                    else:  # Compare synergy of drug actions with interactor/pTAU relations
                        synergy = self.__check_rel_synergy(mapped_actions, rels)
                        if synergy is True:  # Good comparison
                            self.drug_scores[drug_name][INTERACTORS][target_name]['synergy'] = synergy
                            pts += self.__reward

                        else:
                            pts += self.__penalty

                self.drug_scores[drug_name][INTERACTORS][target_name][POINTS] = pts

    @staticmethod
    def __query_db_ct_data(db_id: str) -> Optional[dict]:
        """Checks if query results are stored in cache."""
        logger.info("Querying SQLite DB for Clinical Trial data")

        results = session().query(Trials).filter_by(drugbank_id=db_id).all()
        compiled = dict()
        if results:
            for r in results:
                if r.trial_id is not None:  # No 'trial_id' indicates no associated trials and was previously checked
                    compiled[r.trial_id] = {'trial_status': r.trial_status,
                                            'conditions': r.conditions.split(";"),
                                            'drugs_in_trial': r.drugs_in_trial.split("|")}

            return compiled

    @staticmethod
    def __import_ct_data(ct_data: list):
        """Imports previously missing clinical trial data into SQLite DB."""
        logger.info("Importing missing Clinical Trial data")
        sess = session()
        for entry in ct_data:
            for key, vals in entry.items():
                if isinstance(vals, list):
                    entry[key] = "|".join(vals)

            trial = Trials(**entry)
            sess.add(trial)

        sess.commit()

    def __collect_ct_info(self):
        """Collects clinical trial information for every identified drug."""
        logger.info("Collecting Clinical Trial information for each drug")
        to_import = []
        db_ids = {data[IDENTIFIERS]['drugbank_id']: drug_name for drug_name, data in self.drug_metadata.items()}
        for db_id, drug_name in tqdm(db_ids.items(), total=len(db_ids), desc="Gathering clinical trial info"):
            cached_data = self.__query_db_ct_data(db_id)
            if cached_data is not None:
                self.drug_metadata[drug_name][CLINICAL_TRIALS] = \
                    {**self.drug_metadata[drug_name][CLINICAL_TRIALS], **cached_data}

            else:  # Have to query graphstore
                results = rest_query.sql(CLINICAL_TRIAL_FROM_DRUG.format(db_id)).data
                if results:
                    ct_data = dict()
                    for r in results:
                        # Parse graphstore data
                        status = r['overall_status']
                        trial_id = r['trial_id']
                        primary_condition = r['condition'] if r['condition'] is not None else []
                        mesh = r['mesh_conditions'] if r['mesh_conditions'] is not None else []
                        conditions = ";".join(set(primary_condition + mesh))
                        trial_drugs = r['drugs_in_trial'] if 'drugs_in_trial' in r else None

                        # Structure data
                        import_data = {'trial_id': trial_id, 'drug_name': drug_name, 'drugbank_id': db_id}
                        metadata = {'trial_status': status, 'conditions': conditions, 'drugs_in_trial': trial_drugs}
                        to_import.append({**import_data, **metadata})

                        ct_data[trial_id] = metadata

                    self.drug_metadata[drug_name][CLINICAL_TRIALS] = \
                        {**self.drug_metadata[drug_name][CLINICAL_TRIALS], **ct_data}

                else:  # Enter empty row for drug to indicate no clinical trials associated
                    to_import.append({'drug_name': drug_name, 'drugbank_id': db_id})

        if to_import:
            self.__import_ct_data(to_import)

    def score_cts(self, disease_keyword: str = "Alzheimer Disease", similar_diseases: list = []):
        """Scores drugs based on involvement in a clinical trial. Returned dictionary contains
        points and relevant AD-associated CT information."""

        self.__collect_ct_info()
        logger.info("Scoring Clinical Trial data")
        for drug_name, metadata in self.drug_metadata.items():
            pts = self.__reward  # Default to reward unless changed

            if metadata[CLINICAL_TRIALS]:  # There is clinical trial data
                self.drug_scores[drug_name][CLINICAL_TRIALS]['trials'] = dict()
                for trial_id, trial_data in metadata[CLINICAL_TRIALS].items():
                    ct_score = {'keyword_disease_investigated': False,
                                'trial_ongoing': False,
                                'similar_disease_investigated': False}

                    if disease_keyword in trial_data['conditions']:  # Disease-associated CT
                        ct_score['keyword_disease_investigated'] = True
                        if trial_data['trial_status'] in CT_MAPPER['ongoing']:
                            ct_score['trial_ongoing'] = True
                            pts += self.__penalty

                    else:  # Not associated with primary disease
                        if set(trial_data['conditions']) & set(similar_diseases):  # Trial for similar disease
                            ct_score['similar_disease_investigated'] = True
                            pts += self.__reward * 2  # Double points
                        else:
                            if trial_data['trial_status'] in CT_MAPPER['ongoing']:
                                ct_score['trial_ongoing'] = True
                                pts += self.__penalty

                    self.drug_scores[drug_name][CLINICAL_TRIALS]['trials'][trial_id] = ct_score
            self.drug_scores[drug_name][CLINICAL_TRIALS][POINTS] = pts

    def count_associated_pathways(self) -> dict:
        """Goes through each interactor and checks if it is associated with a KEGG or Pathway Commons pathway."""
        pathways = dict()
        for protein in tqdm(self.interactors, desc="Checking pathways"):
            r = rest_query.sql(ASSOCIATED_PATHWAYS.format(protein, protein)).data
            pathways[protein] = r[0]['count']

        return pathways

    @staticmethod
    def __query_db_edge_counts(symbol: str) -> Optional[dict]:
        """Obtains counts from SQLite DB for given gene symbol."""
        logger.info("Querying SQLite DB for edge counts")
        results = session().query(Edges).filter_by(symbol=symbol).all()
        assert len(results) < 2
        if results:
            hit = results[0]
            data = {'symbol': hit.symbol, 'out_count': hit.out_count,
                    'in_count': hit.in_count, 'both_count': hit.both_count}
            return data

    @staticmethod
    def __query_graphstore_edge_counts(symbol: str) -> dict:
        """Queries the graphstore for the edge counts for given gene symbol."""
        in_count = rest_query.sql(IN_COUNT.format(symbol)).data[0]['number']
        out_count = rest_query.sql(OUT_COUNT.format(symbol)).data[0]['number']
        both_count = out_count + in_count

        counts = {'in_count': in_count, 'out_count': out_count, 'both_count': both_count, 'symbol': symbol}

        logger.info(f"Importing edge counts into DB for {symbol}")
        sess = session()
        edge_entry = Edges(**counts)
        sess.add(edge_entry)
        sess.commit()

        return counts

    def count_edges(self) -> dict:
        """Counts the number of incoming, outgoing, and total edges for each interactor."""
        edge_counts = dict()
        for symbol in tqdm(self.interactor_metadata.keys(), total=len(self.interactor_metadata.keys()),
                           desc="Counting edges"):
            counts = self.__query_db_edge_counts(symbol)
            if not counts:
                counts = self.__query_graphstore_edge_counts(symbol)
            self.interactor_metadata[symbol]['edges'] = counts
            edge_counts[symbol] = counts

        return edge_counts

    @staticmethod
    def __check_bioassay_cache() -> dict:
        """Checks if BioAssay cache file present. If not, creates it. Returns current content."""
        # TODO: Move to SQLite DB
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

    def count_bioassays(self) -> dict:
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

        for symbol, bioassay_count in counts.items():
            self.interactor_metadata[symbol]['bioassays'] = bioassay_count

        return counts
