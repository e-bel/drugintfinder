"""Useful general methods."""
import logging
import pandas as pd

logger = logging.getLogger(__name__)


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def export_table(results: pd.DataFrame, file_path: str):
    """Export results dataframe to path."""
    extension = file_path.split(".")[-1]
    acceptable_file_types = ("csv", "xlsx", "tsv")

    if extension not in acceptable_file_types:
        raise ValueError("Results output path must end with 'xlsx', 'csv', or 'tsv'")

    if extension == "xlsx":
        results.to_excel(file_path, index=False)

    elif extension in ("csv", "tsv"):
        separator = "," if extension == "csv" else "\t"
        results.to_csv(file_path, sep=separator, index=False)

    logger.info(f"Results written to {file_path}")
