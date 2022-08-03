"""Calculate the ranking of the hits."""

import json
import logging
import requests
import pandas as pd

from tqdm import tqdm
from typing import Optional
from ebel_rest import query as rest_query
from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound

from drugintfinder.constants import PATENTS, INTERACTORS, PRODUCTS, IDENTIFIERS, NEGREG, POSREG, SYNERGY, POINTS, \
    DAC, ACTION_MAPPER, CT_MAPPER, CLINICAL_TRIALS, TRIALS, IN_COUNT, OUT_COUNT, PUBCHEM_BIOASSAY_API, UNIPROT_ID, \
    TIC, ASSOCIATED_PATHWAYS, TARGET_COUNT
from drugintfinder.finder import InteractorFinder
from drugintfinder.defaults import session, SIMILAR_DISEASES
from drugintfinder.models import Trials, TargetStats, Patents, Products, Drugs, BioAssays

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Ranker:
    """Rank and annotate druggable interactors identified using the InteractorFinder class."""

    def __init__(
            self,
            name: str,
            node_type: str = "protein",
            pmods: list = None,
            penalty: int = -1,
            reward: int = 1,
            disease_keyword: str = "Alzheimer Disease",
            similar_diseases: list = SIMILAR_DISEASES,
            print_sql: bool = False
    ):
        """Init for Ranker.

        Should be initialized with an InteractorFinder object that has results saved.
        """
        self.__session = session()
        self.node_name = name
        self.pmods = pmods
        self.disease = disease_keyword
        self.similar_diseases = similar_diseases

        self.__finder = InteractorFinder(
            node_name=name, node_type=node_type, pmods=pmods, neighbor_edge_type='causal', print_sql=print_sql
        )
        self.__finder.druggable_interactors()
        self.table = self.__finder.results
        if self.table is None:
            raise ValueError(f"No druggable interactors were found for {self.node_name} + {self.pmods}")

        self.__penalty = penalty
        self.__reward = reward

        self.interactor_metadata = {int_name: dict() for int_name in self.interactors}

        self.__cached_drugs = self.__session.query(Drugs.id,
                                                   Drugs.drug_name,
                                                   Drugs.num_targets,
                                                   Drugs.clinical_trials)\
            .filter(Drugs.drug_name.in_(self.interactor_drugs))\
            .all()

        target_mapper, drug_ids, cts = self.__parse_cached_data()
        self.target_count_mapper = target_mapper
        self.relevant_drug_row_ids = drug_ids
        self.relevant_cts = cts

        self.drug_scores = self.__generate_ranking_dict()
        self.drug_metadata = self.__compile_drug_metadata()

        self.ct_summary = pd.DataFrame()

    def __parse_cached_data(self) -> tuple:
        target_mapper = dict()
        drug_ids = []
        cts = dict()
        for drug_entry in self.__cached_drugs:
            drug_id, drug_name, targets, cts_str = drug_entry
            target_mapper[drug_name] = targets
            drug_ids.append(drug_id)
            cts[drug_name] = [trial_id for trial_id in cts_str.split("|") if trial_id]

        return target_mapper, drug_ids, cts

    @property
    def interactors(self) -> list:
        """Compile a list of all identified interactors."""
        return list(self.__finder.unique_interactors())

    @property
    def interactor_drugs(self) -> list:
        """Compile a list of all identified drugs."""
        return list(self.__finder.unique_drugs())

    @property
    def unique_drug_target_combos(self) -> set:
        """Create a list of tuples uniue combinations of drugs and their targets from the ranking results."""
        unique_pairs = set()
        pairs = self.__finder.drug_and_interactors().to_dict('records')

        for combo in pairs:
            unique_pairs.add((combo['drug'], combo['interactor_name']))

        return unique_pairs

    @property
    def dbid_drugname_mapper(self):
        """Map DrugBank IDs to their drug names."""
        return {data[IDENTIFIERS]['drugbank_id']: drug_name for drug_name, data in self.drug_metadata.items()}

    def rank(self):
        """Score drugs and their targets."""
        self.score_drugs()
        self.score_interactors()

    def score_drugs(self):
        """Parse raw drug metadata and calculate points for each ranking criteria."""
        self.score_drug_relationships()
        self.score_patents_and_products()
        self.score_cts()
        # self.score_homologs()

    def score_interactors(self):
        """Parse raw interactor metadata and calculate points for each ranking criteria."""
        self.count_bioassays()
        self.count_edges()

    def __compile_drug_metadata(self) -> dict:
        metadata = self.__generate_ranking_dict()
        logger.info("Compiling drug metadata")

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
                rel_data = {'relation_type': {r['relation_type']},  # Collect rel types in set for later comparison
                            'actions': r['drug_rel_actions'].split("|") if r['drug_rel_actions'] else None}
                metadata[drug_name][INTERACTORS][interactor_name] = rel_data

        return metadata

    def __generate_ranking_dict(self) -> dict:
        """Compile a generic empty dict to fill with raw and scored data."""
        empty_dict = {drug_name: {PATENTS: dict(),
                                  PRODUCTS: dict(),
                                  IDENTIFIERS: dict(),
                                  INTERACTORS: dict(),
                                  CLINICAL_TRIALS: dict(),
                                  TARGET_COUNT: None,
                                  }
                      for drug_name in self.interactor_drugs}
        return empty_dict

    def score_patents_and_products(self):
        """Parse patent and product information and generates a score."""
        logger.info("Scoring patent and product information")

        relevant_patents = self.__session.query(Drugs.drug_name, Patents.expired) \
            .join(Drugs).filter(Patents.drug_id.in_(self.relevant_drug_row_ids)).all()

        relevant_products = self.__session.query(Drugs.drug_name, Products.has_approved_generic) \
            .join(Drugs).filter(Products.drug_id.in_(self.relevant_drug_row_ids)).all()

        patent_mapper = dict()
        for drug_name, expired in relevant_patents:
            if drug_name in patent_mapper:
                patent_mapper[drug_name].append(expired)
            else:
                patent_mapper[drug_name] = [expired]

        product_mapper = dict()
        for drug_name, approved_generic in relevant_products:
            if drug_name in product_mapper:
                product_mapper[drug_name].append(approved_generic)
            else:
                product_mapper[drug_name] = [approved_generic]

        for drug_name in self.drug_scores.keys():
            all_patents_expired = all(patent_mapper[drug_name]) if drug_name in patent_mapper else False
            approved_generic_available = any(product_mapper[drug_name]) if drug_name in product_mapper else False
            self.drug_scores[drug_name][PATENTS] = {'expired': all_patents_expired,
                                                    POINTS: self.__reward if all_patents_expired else self.__penalty}
            self.drug_scores[drug_name][PRODUCTS] = {'has_approved_generic': approved_generic_available,
                                                     POINTS: self.__reward if approved_generic_available
                                                     else self.__penalty}

            try:
                self.drug_scores[drug_name][TARGET_COUNT] = self.target_count_mapper[drug_name]

            except KeyError:  # If drug not in cache then probably forgot to populate
                logger.error(f"{drug_name} not found in cache. Try populate.populate() again.")

    def score_homologs(self, tc_json: str, db_id_mapper: dict, tc_threshold: int = 0.95):
        """Produce a dictionary detailing structural homologs of every drugbank drug."""
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
        """Check if both 'positive_regulator' and 'negative_regulator' in the set of mapped drug/target relations."""
        contradiction = False
        if POSREG and NEGREG in drug_actions:
            contradiction = True
        return contradiction

    @staticmethod
    def __check_rel_synergy(drug_actions: set, int_rel: set) -> Optional[bool]:
        """Return True is drug/target relation and interactor/pTAU relation results in decrease of pTAU else False."""
        pos_int_rel = ['increases', 'directly_increases']
        neg_int_rel = ['decreases', 'directly_decreases']

        # TODO make generic so it works for any target
        # drug -inhibits-> interactor -decreases-> pTAU: BAD
        if NEGREG in drug_actions and any(rel in int_rel for rel in neg_int_rel):
            return False

        # drug -inhibits-> interactor -increases-> pTAU: GOOD
        elif NEGREG in drug_actions and any(rel in int_rel for rel in pos_int_rel):
            return True

        # drug -activates-> interactor -increases-> pTAU: BAD
        elif POSREG in drug_actions and any(rel in int_rel for rel in pos_int_rel):
            return False

        # drug -activates-> interactor -decreases-> pTAU: GOOD
        elif POSREG in drug_actions and any(rel in int_rel for rel in neg_int_rel):
            return True

        else:  # drug_actions = 'neutral'
            return None

    def score_drug_relationships(self):
        """Add the score to the drug/target pair.

        Score is based on whether it has information on the action (inhibitor/activator) and whether it makes
        sense based on the relationship of the interactor with target.

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
                    SYNERGY: None,
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
                        self.drug_scores[drug_name][INTERACTORS][target_name][SYNERGY] = synergy

                        if synergy is True:  # Good comparison
                            pts += self.__reward

                        else:
                            pts += self.__penalty

                self.drug_scores[drug_name][INTERACTORS][target_name][POINTS] = pts

    @staticmethod
    def __chunk_cts(trial_id_list: list, chunk_size: int = 100) -> list:
        """Create a generator of trial ID lists."""
        total_num_chunks = (len(trial_id_list) // chunk_size) + 1
        chunk_index = 0
        while chunk_index < total_num_chunks:
            start = chunk_index * chunk_size
            stop = (chunk_index + 1) * chunk_size
            yield trial_id_list[start:stop]
            chunk_index += 1

    def __build_ct_cache(self) -> dict:
        """Build a dictionary of metadata for each clinical trial relevant to the query."""
        cached_ct_data = dict()
        relevant_trial_ids = list({trial_id for id_list in self.relevant_cts.values() for trial_id in id_list})

        trial_id_chunks = self.__chunk_cts(relevant_trial_ids, chunk_size=400)
        for id_batch in trial_id_chunks:
            ct_results = self.__session.query(Trials.trial_id, Trials.conditions, Trials.trial_status)\
                .filter(Trials.trial_id.in_(id_batch)).all()

            for ct_entry in ct_results:
                trial_id, conditions, trial_status = ct_entry

                metadata = {'conditions': conditions.split(";"), 'trial_status': trial_status}
                cached_ct_data[trial_id] = metadata

        return cached_ct_data

    def __compile_ct_metadata(self) -> dict:
        """Compile the Clinical Trial data from the DB."""
        ct_mapper = dict()
        cached_ct_data = self.__build_ct_cache()

        for drug_name, trial_ids in self.relevant_cts.items():
            for trial_id in trial_ids:
                metadata = cached_ct_data[trial_id]

                if drug_name in ct_mapper:
                    ct_mapper[drug_name][trial_id] = metadata

                else:
                    ct_mapper[drug_name] = {trial_id: metadata}

        return ct_mapper

    def score_cts(self):
        """Score drugs based on involvement in a clinical trial.

        Returned dictionary contains points and relevant AD-associated CT information.
        """
        logger.info("Scoring Clinical Trial data")
        ct_mapper = self.__compile_ct_metadata()
        ct_summary_rows = []

        for drug_name, ct_metadata in ct_mapper.items():
            pts = self.__reward  # Default to reward unless changed

            self.drug_scores[drug_name][CLINICAL_TRIALS][TRIALS] = dict()
            for trial_id, trial_data in ct_metadata.items():
                ct_score = {'keyword_disease_investigated': False,
                            'trial_ongoing': False,
                            'similar_disease_investigated': False,
                            'conditions_investigated': trial_data['conditions']}

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
                ct_row = {'drug': drug_name, 'trial_id': trial_id, **ct_score}
                ct_summary_rows.append(ct_row)

            self.drug_scores[drug_name][CLINICAL_TRIALS][POINTS] = pts

        # Compile summary dataframe of CT information
        column_names = ["Drug", "Trial ID", "Keyword Disease Investigated", "Trial Ongoing",
                        "Similar Disease Investigated", "Conditions Investigated"]

        if ct_summary_rows:
            self.ct_summary = pd.DataFrame(ct_summary_rows)
            self.ct_summary.columns = column_names
        else:
            self.ct_summary = pd.DataFrame(columns=column_names)

    def count_associated_pathways(self) -> dict:
        """Go through each interactor and checks if it is associated with a KEGG or Pathway Commons pathway."""
        pathways = dict()
        for protein in tqdm(self.interactors, desc="Checking pathways"):
            r = rest_query.sql(ASSOCIATED_PATHWAYS.format(protein, protein)).data
            pathways[protein] = r[0]['count']

        return pathways

    @staticmethod
    def __query_db_edge_counts(symbol: str) -> Optional[dict]:
        """Obtain counts from SQLite DB for given gene symbol."""
        logger.info("Querying SQLite DB for edge counts")
        try:
            results = session().query(TargetStats).filter_by(symbol=symbol).one()

        except MultipleResultsFound as mrfe:
            print(mrfe)
            logger.warning(f"Multiple results found for {symbol}")
            results = None

        except NoResultFound:
            logger.info(f"No edge counts found for {symbol}")
            results = None

        if results:
            data = {'symbol': results.symbol,
                    'out_count': results.out_count,
                    'in_count': results.in_count,
                    'both_count': results.both_count}
            return data

    @staticmethod
    def __query_graphstore_edge_counts(symbol: str) -> dict:
        """Query the graphstore for the edge counts for given gene symbol."""
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
        """Count the number of incoming, outgoing, and total edges for each interactor."""
        edge_counts = dict()
        logger.info("Gathering edge counts")
        for symbol in self.interactor_metadata.keys():
            counts = self.__query_db_edge_counts(symbol)
            if not counts:
                counts = self.__query_graphstore_edge_counts(symbol)
            self.interactor_metadata[symbol]['edges'] = counts
            edge_counts[symbol] = counts

        return edge_counts

    def __get_bioassay_cache(self) -> dict:
        """Get BioAssay counts in DB and returns a dict - key: protein target symbol, val: bioassay count."""
        bioassay_counts = self.__session.query(BioAssays.target, BioAssays.num_assays).all()
        return {r[0]: r[1] for r in bioassay_counts}

    def count_bioassays(self) -> dict:
        """Query the PubChem API to obtain the number of available BioAssays for each gene symbol in the list.

        Returns
        -------
        dict
            Key is the gene symbol, value is the number of BioAssays available.
        """
        counts = self.__get_bioassay_cache()
        symbols_missing_data = set(self.interactors) - set(counts.keys())  # Create set of symbols not in cache

        if symbols_missing_data:
            for symbol in tqdm(symbols_missing_data, desc="Counting BioAssays for targets"):

                up_acc_results = rest_query.sql(UNIPROT_ID.format(symbol)).data

                if up_acc_results:
                    up_acc = up_acc_results[0]['uniprot.id']

                    filled = PUBCHEM_BIOASSAY_API.format(up_acc)
                    resp = requests.get(filled)
                    number_assays = len(resp.text.split("\n")) - 1  # Remove header
                    counts[symbol] = number_assays

                    row = BioAssays(target=symbol, num_assays=number_assays)
                    self.__session.add(row)

            self.__session.commit()

        for symbol in self.interactor_metadata.keys():
            self.interactor_metadata[symbol]['bioassays'] = counts[symbol]  # Number of assays

        return counts

    def __create_summary_table(self) -> pd.DataFrame:
        rows = []
        for drug_name, symbol in self.unique_drug_target_combos:
            drug_entry = self.drug_scores[drug_name]

            synergy_result = drug_entry[INTERACTORS][symbol][SYNERGY]
            if synergy_result is None:
                synergizes = "N/A"

            elif synergy_result is False:
                synergizes = "No"

            else:  # is True
                synergizes = "Yes"

            # Interactor metadata
            num_bioassays = self.interactor_metadata[symbol]['bioassays']
            num_total_edges = self.interactor_metadata[symbol]['edges']['both_count']

            # Drug metadata
            ongoing_patent = "Yes" if drug_entry[PATENTS]['expired'] is True else "No"
            has_generic = "Yes" if drug_entry[PRODUCTS]['has_approved_generic'] is True else "No"
            target_count = drug_entry[TARGET_COUNT]
            target_count_entry = "N/A" if target_count == -1 or not target_count else target_count

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

    def summarize(self, pivot: bool = False) -> pd.DataFrame:
        """Summarize the ranking results.

        Be default, the summarized table is compiled into rows of unique drug/targets. If `pivot` is enabled, the
        table instead focuses on the targets and the columns adjust accordingly.

        :param pivot: If True, drugs are compressed into a new column and table now focuses on targets.
        :return: pd.DataFrame os summary results.
        """
        summary_df = self.__create_summary_table()

        if pivot:
            num_drugs_mapper = summary_df.Target.value_counts()
            summary_df = summary_df.drop([
                'Drug',
                'Synergizes',
                'Drug Patent Ongoing',
                'Generic Version of Drug Available',
                'Number of Drug Targets'
            ], axis=1).drop_duplicates().reset_index(drop=True)
            summary_df["Number of Known Drugs"] = [num_drugs_mapper[target] for target in summary_df.Target]
            summary_df["BioAssay to Known Drug Ratio"] = \
                summary_df["Number of BioAssays for Target"] / summary_df["Number of Known Drugs"]
            summary_df = summary_df.sort_values(by="Number of BioAssays for Target", ascending=False)

        return summary_df
