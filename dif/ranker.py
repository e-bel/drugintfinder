"""Calculate the ranking of the hits."""

import os
import json
import requests
import logging
import pandas as pd

from tqdm import tqdm
from typing import Optional
from datetime import datetime
from json.decoder import JSONDecodeError
from ebel_rest import query as rest_query

from dif.constants import *
from dif.finder import InteractorFinder
from dif.defaults import BIOASSAY_CACHE, session, SIMILAR_DISEASES
from dif.models import Trials, TargetStats, Patents, Products, Drugs

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Ranker:
    """Class for taking an InteractorFinder object with saved results and annotating said results with ranking
    metadata."""

    def __init__(self, symbol: str, pmods: list = None, penalty: int = -1, reward: int = 1,
                 disease_keyword: str = "Alzheimer Disease", similar_diseases: list = SIMILAR_DISEASES):
        """Should be initialized with an InteractorFinder object that has results saved."""
        self.__session = session()
        self.symbol = symbol
        self.pmods = pmods
        self.disease = disease_keyword
        self.similar_diseases = similar_diseases

        self.__finder = InteractorFinder(symbol=symbol, pmods=pmods, edge='causal')
        self.__finder.druggable_interactors()
        self.table = self.__finder.results
        if self.table is None:
            raise ValueError(f"No druggable interactors were found for {self.symbol} + {self.pmods}")

        self.__penalty = penalty
        self.__reward = reward

        self.interactor_metadata = {int_name: dict() for int_name in self.interactors}

        self.drug_metadata = self.__compile_drug_metadata()
        self.drug_scores = self.__generate_ranking_dict()

    @property
    def drug_cache(self) -> dict:
        cached_drugs = self.__session.query(Drugs.drug_name, Drugs.id, Drugs.num_targets).all()
        return {x[0]: {'row_id': x[1], 'num_targets': x[2]} for x in cached_drugs}

    @property
    def relevant_drug_row_ids(self) -> list:
        return [self.drug_cache[drug_name] for drug_name in self.interactor_drugs]

    @property
    def interactors(self) -> list:
        return self.__finder.unique_interactors()

    @property
    def interactor_drugs(self) -> list:
        return self.__finder.unique_drugs()

    @property
    def unique_drug_target_combos(self) -> set:
        unique_pairs = set()
        pairs = self.__finder.drug_and_interactors().to_dict('records')

        for combo in pairs:
            unique_pairs.add((combo['drug'], combo['interactor_name']))

        return unique_pairs

    @property
    def dbid_drugname_mapper(self):
        return {data[IDENTIFIERS]['drugbank_id']: drug_name for drug_name, data in self.drug_metadata.items()}

    def rank(self):
        self.score_drugs()
        self.score_interactors()

    def score_drugs(self):
        """Wrapper method to parse raw drug metadata and calculate points for each ranking criteria."""
        self.score_drug_relationships()
        self.score_ppt()
        self.score_cts()
        # self.score_homologs()

    def score_interactors(self):
        """Wrapper method to parse raw interactor metadata and calculate points for each ranking criteria."""
        self.count_bioassays()
        self.count_edges()

    def __query_db_pp_drugs(self) -> tuple:
        """Gets currently cached patent and product information and adds points."""
        patents, products, targets = dict(), dict(), dict()
        patent_rows = session().query(Patents).all()
        product_rows = session().query(Products).all()

        for pat in patent_rows:
            patents[pat.drug_name] = {'has_patent': pat.has_patent,
                                      'expired': pat.expired,
                                      'patent_numbers': pat.patent_numbers.split("|"),
                                      POINTS: self.__reward if pat.expired else self.__penalty}

        for prod in product_rows:
            products[prod.drug_name] = {'has_generic': prod.has_generic,
                                        'has_approved_generic': prod.has_approved_generic,
                                        'generic_products': prod.generic_products.split("|") or None,
                                        POINTS: self.__reward if prod.has_approved_generic else self.__penalty}

        return patents, products, targets

    def score_ppt(self):
        """Scores the patents, products, and counts targets for each drug."""
        missing_pp_db_ids = set()
        patents, products, targets = self.__query_db_pp_drugs()

        for db_id, drug_name in self.dbid_drugname_mapper.items():
            # Check patents/products - if present then add directly to scores
            if drug_name in patents:
                self.drug_scores[drug_name][PATENTS] = patents[drug_name]

            if drug_name in products:
                self.drug_scores[drug_name][PRODUCTS] = products[drug_name]

            if drug_name in targets:
                self.drug_scores[drug_name][TARGET_COUNT] = targets[drug_name]

            # If one is missing, add drugbank ID to list to be queried in graphstore
            if drug_name not in patents or drug_name not in products or drug_name not in targets:
                missing_pp_db_ids.add(db_id)

        # Query graphstore and parse results
        self.score_patents()
        self.score_products()

    def __compile_patent_information(self) -> dict:
        """Gathers all requested patent information into drug_metadata."""
        relevant_patents = self.session.query(Drugs.drug_name, Patents.expired)\
            .join(Drugs).filter(Patents.drug_id.in_(self.relevant_drug_row_ids)).all()

        patent_mapper = dict()
        for drug_name, expired in relevant_patents:
            if drug_name in patent_mapper:
                patent_mapper[drug_name].append(expired)
            else:
                patent_mapper[drug_name] = [expired]

        return {drug_name: all(expired) for drug_name, expired in patent_mapper}

    def __compile_product_information(self) -> dict:
        """Gathers all requested product information into drug_metadata."""
        relevant_products = self.session.query(Drugs.drug_name, Products.has_approved_generic)\
            .join(Drugs).filter(Products.drug_id.in_(self.relevant_drug_row_ids)).all()

        product_mapper = dict()
        for drug_name, approved_generic in relevant_products:
            if drug_name in product_mapper:
                product_mapper[drug_name].append(approved_generic)
            else:
                product_mapper[drug_name] = [approved_generic]

        return {drug_name: all(approved_generic) for drug_name, approved_generic in product_mapper}

    def __compile_drug_metadata(self) -> dict:
        metadata = self.__generate_ranking_dict()

        for r in self.table.to_dict('records'):
            drug_name = r['drug']
            db_id = r['drugbank_id']
            interactor_name = r['interactor_name']
            metadata[drug_name][IDENTIFIERS] = {'drugbank_id': db_id,
                                                'chembl_id': r['chembl_id'],
                                                'pubchem_id': r['pubchem_id']}

            # Parse relations
            if r['relation_type'] == 'regulates':
                continue

            elif interactor_name in metadata[drug_name][INTERACTORS]:
                metadata[drug_name][INTERACTORS][interactor_name]['relation_type'].add(r['relation_type'])

            else:
                rel_data = {'relation_type': {r['relation_type']},  # Collect relation types in set for later comparison
                            'actions': r['drug_rel_actions'].split("|") if r['drug_rel_actions'] else None}
                metadata[drug_name][INTERACTORS][interactor_name] = rel_data

            metadata[drug_name][TARGET_COUNT] = self.drug_cache[drug_name]['num_targets']

            patent_metadata = self.__compile_patent_information()
            product_metadata =

        return metadata

    def __generate_ranking_dict(self) -> dict:
        """Compiles a generic empty dict to fill with raw and scored data."""
        empty_dict = {drug_name: {PATENTS: dict(),
                                  PRODUCTS: dict(),
                                  IDENTIFIERS: dict(),
                                  INTERACTORS: dict(),
                                  CLINICAL_TRIALS: dict(),
                                  TARGET_COUNT: None,
                                  }
                      for drug_name in self.interactor_drugs}
        return empty_dict

    def score_patents(self):
        """Parses patent information from graphstore and generates a score. Imports into SQLite DB at the end."""
        drugs_and_patents_raw = {drug_name: dd[PATENTS] for drug_name, dd in self.drug_metadata.items()}
        logger.info("Scoring patent information")
        for drug_name, patent_info in drugs_and_patents_raw.items():
            patent_numbers = []
            expired_list = []
            if not self.drug_scores[drug_name][PATENTS]:  # No cached data

                    all_expired = all(expired_list)
                    import_data = {'has_patent': bool(patent_numbers),
                                   'expired': all_expired,
                                   'patent_numbers': "|".join(patent_numbers),
                                   POINTS: self.__reward if all_expired else self.__penalty}

                else:
                    import_data = {'has_patent': False, 'expired': False,
                                   'patent_numbers': "N/A", POINTS: self.__penalty}

                self.drug_scores[drug_name][PATENTS] = import_data
        # TODO finish

    def score_products(self):
        """Gives scores to drugs in the hit list based on whether they have an approved generic version."""
        drugs_and_product_info = {drug_name: dd[PRODUCTS] for drug_name, dd in self.drug_metadata.items()}
        logger.info("Scoring product information")
        for drug_name, product_info in drugs_and_product_info.items():
            if not self.drug_scores[drug_name][PRODUCTS]:  # No cached data
                if has_generic and has_approved_generic:
                    pts = self.__reward
                else:
                    pts = self.__penalty
        # TODO finish

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
        logger.info("Scoring drug relationships")
        for drug_name, metadata in self.drug_metadata.items():
            for target_name, ti_metadata in metadata[INTERACTORS].items():
                self.drug_scores[drug_name][INTERACTORS][target_name] = {
                    TIC: False,
                    DAC: False,
                    SYNERGY: False,
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
                            self.drug_scores[drug_name][INTERACTORS][target_name][SYNERGY] = synergy
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

    def score_cts(self):
        """Scores drugs based on involvement in a clinical trial. Returned dictionary contains
        points and relevant AD-associated CT information."""

        self.__collect_ct_info()
        logger.info("Scoring Clinical Trial data")
        for drug_name, metadata in self.drug_metadata.items():
            pts = self.__reward  # Default to reward unless changed

            if metadata[CLINICAL_TRIALS]:  # There is clinical trial data
                self.drug_scores[drug_name][CLINICAL_TRIALS][TRIALS] = dict()
                for trial_id, trial_data in metadata[CLINICAL_TRIALS].items():
                    ct_score = {'keyword_disease_investigated': False,
                                'trial_ongoing': False,
                                'similar_disease_investigated': False}

                    if self.disease in trial_data['conditions']:  # Disease-associated CT
                        ct_score['keyword_disease_investigated'] = True
                        if trial_data['trial_status'] in CT_MAPPER['ongoing']:
                            ct_score['trial_ongoing'] = True
                            pts += self.__penalty

                    else:  # Not associated with primary disease
                        if set(trial_data['conditions']) & set(self.similar_diseases):  # Trial for similar disease
                            ct_score['similar_disease_investigated'] = True
                            pts += self.__reward * 2  # Double points
                        else:
                            if trial_data['trial_status'] in CT_MAPPER['ongoing']:
                                ct_score['trial_ongoing'] = True
                                pts += self.__penalty

                    self.drug_scores[drug_name][CLINICAL_TRIALS][TRIALS][trial_id] = ct_score
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
        results = session().query(TargetStats).filter_by(symbol=symbol).all()
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
        edge_entry = TargetStats(**counts)
        sess.add(edge_entry)
        sess.commit()

        return counts

    def count_edges(self) -> dict:
        """Counts the number of incoming, outgoing, and total edges for each interactor."""
        edge_counts = dict()
        logger.info("Gathering edge counts")
        for symbol in tqdm(self.interactor_metadata.keys(), total=len(self.interactor_metadata),
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

    def summarize(self) -> pd.DataFrame:
        rows = []
        for drug_name, symbol in self.unique_drug_target_combos:

            synergizes = "Yes" if self.drug_scores[drug_name][INTERACTORS][symbol][SYNERGY] else "No"

            # Interactor metadata
            num_bioassays = self.interactor_metadata[symbol]['bioassays']
            num_total_edges = self.interactor_metadata[symbol]['edges']['both_count']

            # Drug metadata
            ongoing_patent = "Yes" if self.drug_scores[drug_name][PATENTS]['expired'] is True else "No"
            has_generic = "Yes" if self.drug_scores[drug_name][PRODUCTS]['has_generic'] is True else "No"
            target_count = self.drug_scores[drug_name][TARGET_COUNT]
            target_count_entry = "N/A" if target_count == -1 else target_count

            # Compile
            row = {'Drug': drug_name,
                   'Target': symbol,
                   'Synergizes': synergizes,
                   'Number of BioAssays for Target': num_bioassays,
                   'Number of Causal Edges for Target': num_total_edges,
                   'Drug Patent Ongoing': ongoing_patent,
                   'Generic Version of Drug Available': has_generic,
                   'Number of Drug Targets': target_count_entry}
            rows.append(row)

        return pd.DataFrame(rows)

    def count_number_off_targets(self):
        """Calculates the number of "off targets" per compound i.e. how many compounds it targets."""
        # TODO
        pass
