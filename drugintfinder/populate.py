"""Populate the database with the base information for faster querying."""
import logging

from tqdm import tqdm
from sqlalchemy import func
from datetime import datetime
from typing import List, Union, Dict
from ebel_rest import query as rest_query

from drugintfinder.constants import DRUG_COUNT, DRUG_METADATA, CLINICAL_TRIALS_COUNT, CLINICAL_TRIALS_DATA
from drugintfinder.defaults import session
from drugintfinder.models import Trials, Patents, Products, Drugs

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

sess = session()


def populate():
    """Populate tables."""
    ClinicalTrialPopulator().populate()
    DrugPopulator().populate()


class DrugPopulator:
    """Populate the associated drug tables."""

    def __init__(self):
        """Init method for DrugPopulator."""
        self.__num_drugs_in_graphstore = rest_query.sql(DRUG_COUNT).data[0]['num_drugs']

    def __collect_drugs_from_graphstore(self, chunk_size: int = 10000) -> list:
        """Query clinical trial data in chunks."""
        logger.warning("Collecting drug metadata from graphstore")
        total_num_chunks = (self.__num_drugs_in_graphstore // chunk_size) + 1

        chunk_index = 0
        while chunk_index < total_num_chunks:
            table_chunk = rest_query.sql(DRUG_METADATA.format(chunk_size * chunk_index, chunk_size)).data
            yield table_chunk
            chunk_index += 1

    @staticmethod
    def __extract_values(value_dict: dict, key: str):
        """Extract non-None values from key if present and pops key."""
        values = []
        if key in value_dict:
            vals = value_dict.pop(key)
            if vals is not None:
                for val in vals:
                    if val is not None:
                        values.append(val)

        return value_dict, values

    def populate(self):
        """Populate the SQLite DB with Drug metadata."""
        update = self.__update_needed()
        if update:
            logger.info("Populating drug information")
            self.__parse__and_import_drug_metadata()

    def __update_needed(self) -> bool:
        """Check if drug data is missing."""
        db_entries = sess.query(func.count(Drugs.id)).first()[0]
        return False if db_entries == self.__num_drugs_in_graphstore else True

    def __parse__and_import_drug_metadata(self):
        drug_generator = self.__collect_drugs_from_graphstore()
        logger.info("Parsing and importing drug information from graphstore")
        for drug_data_chunk in drug_generator:
            for drug_entry in tqdm(drug_data_chunk, desc="Parsing and importing drug data"):
                patents = drug_entry.pop("drug_patents")
                products = drug_entry.pop("drug_products")
                drug_entry, targets = self.__extract_values(drug_entry, "target_symbols")

                drug_entry, trial_table_values = self.__extract_values(drug_entry, "clinical_trials")

                patent_rows = self.__parse_patents(patents) if patents else []
                product_rows = self.__parse_products(products) if products else []

                drug_entry['num_targets'] = len(targets)
                drug_entry['targets'] = "|".join(targets)
                drug_entry['clinical_trials'] = "|".join(trial_table_values)

                new_drug = Drugs(**drug_entry)
                new_drug.patents = patent_rows
                new_drug.products = product_rows
                sess.add(new_drug)

        sess.commit()

    @staticmethod
    def __get_clinical_trial_ids() -> dict:
        """Get row entries and Clinical Trial IDs for each CT row in DB."""
        trial_rows = sess.query(Trials, Trials.trial_id).all()
        return {r[1]: r[0] for r in trial_rows}

    @staticmethod
    def __parse_patents(patents: Dict[str, dict]) -> list:
        """Parse patent information from graphstore and imports it into SQLite DB at the end."""
        patent_rows = []
        patent_info = patents['patent']

        if not isinstance(patent_info, list):
            patent_info = [patent_info]

        for indiv_patent in patent_info:

            expired = datetime.today() > datetime.strptime(indiv_patent['expires'], "%Y-%m-%d")  # Boolean
            pat_number = indiv_patent['number']

            import_data = {'expired': expired, 'patent_number': pat_number}
            new_patent = Patents(**import_data)
            patent_rows.append(new_patent)
            sess.add(new_patent)

        return patent_rows

    @staticmethod
    def __parse_products(products_raw: Union[Dict[str, dict], List[Dict[str, dict]]]) -> list:
        """Parse product information from graphstore and imports it into SQLite DB at the end."""
        product_rows = []

        if not isinstance(products_raw, list):
            products_raw = [products_raw]

        labelled_products = {product['name']: {'has_generic': False,
                                               'is_approved': False,
                                               'has_approved_generic': False}
                             for product in products_raw}

        for product in products_raw:  # Usually several disrtibutors for the same drug
            product_name = product['name']

            if product['approved'] == 'true':
                labelled_products[product_name]['is_approved'] = True

            if product['generic'] == 'true':
                labelled_products[product_name]['has_generic'] = True

            if labelled_products[product_name]['is_approved'] and labelled_products[product_name]['has_generic']:
                labelled_products[product_name]['has_approved_generic'] = True

        for product_name, metadata in labelled_products.items():
            import_metadata = {'product_name': product_name, **metadata}

            new_product = Products(**import_metadata)
            product_rows.append(new_product)
            sess.add(new_product)

        return product_rows


class ClinicalTrialPopulator:
    """Clinical Trials class for importing into DB."""

    def __init__(self):
        """Init method to CTPopulator."""
        self.__num_trials_in_graphstore = rest_query.sql(CLINICAL_TRIALS_COUNT).data[0]['trial_count']

    def populate(self):
        """Populate the SQLite DB with ClinicalTrial data."""
        update = self.__update_needed()
        if update:
            logger.info("Populating Clinical Trial information")
            ct_data = self.__collect_ct_info()
            self.__populate_table(ct_data)

    def __update_needed(self) -> bool:
        """Check if Clinical Trial data is missing."""
        db_entries = sess.query(func.count(Trials.id)).first()[0]
        return False if db_entries == self.__num_trials_in_graphstore else True

    def __collect_ct_info_from_graphstore(self, chunk_size: int = 10000) -> list:
        """Query clinical trial data in chunks."""
        total_num_chunks = (self.__num_trials_in_graphstore // chunk_size) + 1

        chunk_index = 0
        while chunk_index < total_num_chunks:
            table_chunk = rest_query.sql(CLINICAL_TRIALS_DATA.format(chunk_size * chunk_index, chunk_size)).data
            yield table_chunk
            chunk_index += 1

    def __collect_ct_info(self) -> list:
        """Collect clinical trial information for every identified drug."""
        num_chunks = (self.__num_trials_in_graphstore // 10000) + 1
        data_generator = self.__collect_ct_info_from_graphstore()

        to_import = []
        for data_chunk in tqdm(data_generator, total=num_chunks, desc="Collecting all CT data"):
            if data_chunk:
                for r in data_chunk:
                    # Parse graphstore data
                    status = r['overall_status']
                    trial_id = r['trial_id']
                    primary_condition = r['condition'] if r['condition'] is not None else []
                    mesh = r['mesh_conditions'] if r['mesh_conditions'] is not None else []
                    conditions = ";".join(set(primary_condition + mesh))
                    trial_drugs = r['drugs_in_trial'] if 'drugs_in_trial' in r else None

                    # Structure data
                    metadata = {'trial_id': trial_id, 'trial_status': status,
                                'conditions': conditions, 'drugs_in_trial': trial_drugs}

                    to_import.append(metadata)

        return to_import

    @staticmethod
    def __populate_table(data: list):
        """Populate the SQLite DB with ClinicalTrial data."""
        logger.info("Importing Clinical Trial data")
        for entry in tqdm(data, desc="Parsing and importing clinical trials data"):
            for key, vals in entry.items():
                if isinstance(vals, list):
                    entry[key] = "|".join(vals)

            row = Trials(**entry)
            sess.add(row)

        sess.commit()
