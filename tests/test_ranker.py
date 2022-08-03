"""Tests for `Ranker` class."""
import pandas as pd
import pytest

from drugintfinder.ranker import Ranker
from .constants import MAPT, PHOSPHORYLATION


@pytest.fixture
def ranker():
    """Create the base Ranker object for testing."""
    ranker = Ranker(name=MAPT, pmods=[PHOSPHORYLATION], print_sql=False)
    return ranker


class TestRanker:
    """Tests for the Ranker class."""

    def test_rank(self, ranker):
        """Test the find_interactors method."""
        ranker.rank()
        summary = ranker.summarize()

        assert isinstance(summary, pd.DataFrame)
        assert len(summary) >= 847

        expected_cols = ['Drug', 'Target', 'Synergizes', 'Number of BioAssays for Target',
                         'Number of Causal Edges for Target', 'Drug Patent Ongoing',
                         'Generic Version of Drug Available', 'Number of Drug Targets']
        assert all([col in summary.columns for col in expected_cols])

# TODO finish writing tests
