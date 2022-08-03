"""Tests for `InteractorFinder` class."""
import pandas as pd
import pytest

from drugintfinder.finder import InteractorFinder
from .constants import CD33, PROTEIN, PHOSPHORYLATION, CAUSAL


@pytest.fixture
def finder():
    """Create reusable interactor finder object."""
    finder = InteractorFinder(
        node_name=CD33, pmods=[PHOSPHORYLATION], neighbor_edge_type=CAUSAL, node_type=PROTEIN, print_sql=False
    )
    return finder


class TestInteractorFinder:
    """Tests for the InteractorFinder class."""

    def test_find_interactors(self, finder):
        """Test the find_interactors method."""
        finder.find_interactors()
        results = finder.results

        assert results is not None
        assert isinstance(results, pd.DataFrame)
        assert len(results) >= 11

        expected_cols = ["target_species", "pmid", "pmc", "interactor_type", "interactor_name", "interactor_bel",
                         "relation_type", "target_bel", "target_type", "target_symbol", "pmod_type"]

        assert all([col in results.columns for col in expected_cols])

    def test_druggable_interactors(self, finder):
        """Test the druggable_interactors method."""
        finder.druggable_interactors()
        results = finder.results

        assert results is not None
        assert isinstance(results, pd.DataFrame)
        assert len(results) >= 70

        expected_cols = ['drug', 'capsule_interactor_type', 'capsule_interactor_bel', 'interactor_bel',
                         'interactor_type', 'interactor_name', 'relation_type', 'target_bel', 'target_symbol',
                         'target_type', 'pmid', 'pmc', 'rel_pub_year', 'rel_rid', 'drug_rel_rid', 'drug_rel_actions',
                         'drugbank_id', 'chembl_id', 'pubchem_id', 'pmod_type']

        assert all([col in results.columns for col in expected_cols])

    def test_unique_interactors(self, finder):
        """Test the unique_interactors method."""
        finder.druggable_interactors()

        ui = finder.unique_interactors()
        assert isinstance(ui, tuple)
        assert len(ui) >= 4

    def test_unique_drugs(self, finder):
        """Test the unique_drugs method."""
        finder.druggable_interactors()

        ud = finder.unique_drugs()
        assert isinstance(ud, tuple)
        assert len(ud) >= 20
