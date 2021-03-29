"""Populate the database with the base information for faster querying."""
import logging

from tqdm import tqdm
from datetime import datetime
from ebel_rest import query as rest_query

from dif.constants import *
from dif.defaults import session
from dif.models import Trials, Patents, Products, Drugs

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def populate():
    populate_ct()
    populate_drugs()


def populate_drugs():
    """Populates the SQLite DB with Drug data."""
    pass


# Drugs
def __collect_drugs_from_graphstore() -> list:
    """Query clinical trial data in chunks"""
    logger.info("Collecting drug metadata from graphstore")
    drug_data = rest_query.sql(DRUG_METADATA).data
    return drug_data


def __collect_drug_info():
    sess = session()
    drug_data = __collect_drugs_from_graphstore()
    trial_mapper = __get_clinical_trial_ids()
    for drug_entry in drug_data:
        patents = drug_entry.pop("drug_patents")
        products = drug_entry.pop("drug_products")
        targets = drug_entry.pop("target_symbols")
        trial_table_ids = [trial_mapper[trial_id] for trial_id in drug_entry.pop("clinical_trials")]

        patent_rows = score_patents(patents)
        product_rows = score_products(products)

        drug_entry['num_targets'] = len(targets)
        drug_entry['targets'] = "|".join(targets)
        drug_entry['clinical_trials'] = trial_table_ids

        new_drug = Drugs(**drug_entry)
        new_drug.patents = patent_rows
        new_drug.products = product_rows
        sess.add(new_drug)

    sess.commit()


def __get_clinical_trial_ids() -> dict:
    """Get row IDs and Clinical Trial IDs for each CT row in DB."""
    trial_rows = session().query(Trials.id, Trials.trial_id).all()
    return {r[1]: r[0] for r in trial_rows}


def score_patents(patents: list) -> list:
    """Parses patent information from graphstore and imports it into SQLite DB at the end."""
    patent_rows = []
    sess = session()
    for patent_info in patents:
        patent_list = patent_info['patent']

        if not isinstance(patent_info['patent'], list):
            patent_list = [patent_info['patent']]

        for indiv_patent in patent_list:
            expired = datetime.today() > datetime.strptime(indiv_patent['expires'], "%Y-%m-%d")  # Boolean
            pat_number = indiv_patent['number']

            import_data = {'expired': expired, 'patent_number': pat_number}
            new_patent = Patents(**import_data)
            patent_rows.append(new_patent)
            sess.add(new_patent)

    sess.commit()
    return patent_rows


def score_products(products_raw: list) -> list:
    """Parses product information from graphstore and imports it into SQLite DB at the end."""
    sess = session()
    product_rows = []
    for product_info in products_raw:
        has_generic = False
        is_approved = False
        has_approved_generic = False
        product_name = None

        if not isinstance(product_info, list):
            product_info = [product_info]

        for product in product_info:
            if product['generic'] == 'true':
                has_generic = True
                is_approved = True if product['approved'] == 'true' else False
                if is_approved:
                    has_approved_generic = True
                product_name = product['name']

        import_metadata = {'has_generic': has_generic,
                           'has_approved_generic': has_approved_generic,
                           'product_name': product_name,
                           'is_approved': is_approved}

        new_product = Products(**import_metadata)
        product_rows.append(new_product)
        sess.add(new_product)

    sess.commit()
    return product_rows


# Clinical Trials
def populate_ct():
    """Populates the SQLite DB with ClinicalTrial data."""
    ct_data = __collect_ct_info()
    __populate_table(Trials, data=ct_data, desc="Clinical Trials")


def __collect_ct_info_from_graphstore(num_trials: int, chunk_size: int = 10000) -> list:
    """Query clinical trial data in chunks"""
    total_num_chunks = (num_trials // chunk_size) + 1

    chunk_index = 0
    while chunk_index < total_num_chunks:
        table_chunk = rest_query.sql(CLINICAL_TRIALS_DATA.format(chunk_size * chunk_index, chunk_size)).data
        yield table_chunk
        chunk_index += 1


def __collect_ct_info() -> list:
    """Collects clinical trial information for every identified drug."""
    logger.info("Collecting Clinical Trial information")
    num_trials = rest_query.sql(CLINICAL_TRIALS_COUNT).data[0]['trial_count']
    num_chunks = (num_trials // 10000) + 1
    data_generator = __collect_ct_info_from_graphstore(num_trials)

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

    return to_import


def __populate_table(model, data: list, desc: str = ""):
    """Populates the SQLite DB with ClinicalTrial data."""
    logger.info(f"Importing {desc} data")
    sess = session()
    for entry in data:
        for key, vals in entry.items():
            if isinstance(vals, list):
                entry[key] = "|".join(vals)

        row = model(**entry)
        sess.add(row)

    sess.commit()
