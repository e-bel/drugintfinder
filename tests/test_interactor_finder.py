"""Tests for `InteractorFinder` class."""
import pandas as pd

from drugintfinder.finder import InteractorFinder
from .constants import MAPT, PROTEIN, PHOSPHORYLATION, CAUSAL


finder = InteractorFinder(symbol=MAPT, pmods=[PHOSPHORYLATION], edge=CAUSAL)


class TestInteractorFinder:
    """Tests for the InteractorFinder class."""

    def test_find_interactors(self):
        """Test the find_interactors method."""
        finder.find_interactors(target_type=PROTEIN)
        results = finder.results

        assert results is not None
        assert isinstance(results, pd.DataFrame)
        assert len(results) > 1500

        expected_cols = ["target_species", "pmid", "pmc", "interactor_type", "interactor_name", "interactor_bel",
                         "relation_type", "target_bel", "target_type", "target_symbol", "pmod_type"]

        assert all([col in results.columns for col in expected_cols])

    def test_druggable_interactors(self):
        """Test the druggable_interactors method."""
        finder.druggable_interactors()
        results = finder.results

        assert results is not None
        assert isinstance(results, pd.DataFrame)
        assert len(results) > 57000

        expected_cols = ['drug', 'capsule_interactor_type', 'capsule_interactor_bel', 'interactor_bel',
                         'interactor_type', 'interactor_name', 'relation_type', 'target_bel', 'target_symbol',
                         'target_type', 'pmid', 'pmc', 'rel_pub_year', 'rel_rid', 'drug_rel_rid', 'drug_rel_actions',
                         'drugbank_id', 'chembl_id', 'pubchem_id', 'pmod_type']

        assert all([col in results.columns for col in expected_cols])

    def test_unique_interactors(self):
        """Test the unique_interactors method."""
        finder.druggable_interactors()

        ui = finder.unique_interactors()
        assert isinstance(ui, tuple)
        assert len(ui) == 80

    def test_unique_drugs(self):
        """Test the unique_drugs method."""
        finder.druggable_interactors()

        ud = finder.unique_drugs()
        assert isinstance(ud, tuple)
        assert len(ud) == 670
