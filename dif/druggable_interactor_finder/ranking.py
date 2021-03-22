"""Calculate the ranking of the hits."""

import os
import requests
import pandas as pd

from tqdm.notebook import tqdm

from dif.constants import PUBCHEM_BIOASSAY_API


def count_bioassays(di_targets: list) -> list:
    """Queries the PubChem API to obtain the number of available BioAssays for each gene symbol in the list.

    Parameters
    ----------
    di_targets: list
        List of druggable interactor gene symbols.

    Returns
    -------
    dict
        Key is the gene symbol, value is the number of BioAssays available.
    """
    counts = []
    for symbol in tqdm(di_targets, desc="Counting BioAssays"):
        up_acc_results = b.query_class(class_name='protein', columns=['uniprot.id'], pure=True, name=symbol)

        if 'uniprot' not in up_acc_results[0]:
            up_acc = up_acc_results[1]['uniprot']

        else:
            up_acc = up_acc_results[0]['uniprot']

        filled = PUBCHEM_BIOASSAY_API.format(up_acc)
        resp = requests.get(filled)
        number_assays = len(resp.text.split("\n")) - 1  # Remove header
        counts.append((symbol, number_assays))
