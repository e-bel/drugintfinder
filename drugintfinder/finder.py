"""Methods for executing pipeline."""
import logging

import pandas as pd

from typing import Optional
from ebel_rest import query as rest_query

from drugintfinder.populate import populate
from drugintfinder.defaults import session
from drugintfinder.models import General, Druggable
from drugintfinder.constants import INTERACTOR_QUERY, PURE_DRUGGABLE_QUERY, CAPSULE_DRUGGABLE_COMPLEX, \
    CAPSULE_DRUGGABLE_MODIFIED, EDGE_MAPPER

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class InteractorFinder:
    """Find the different types of interactors for a given node symbol and protein specification."""

    def __init__(
            self,
            node_name: str,
            node_type: str = 'protein',
            pmods: list = None,
            neighbor_edge_type: str = 'E',
            neighbor_edge_filters: list = None,
            print_sql: bool = False
    ):
        """Init method for InteractorFinder class.

        Parameters
        ----------
        node_name : str
            Value of target node for which you want to find interactors for. This could be a gene symbol, disease name,
            or any other accepted BEL value.
        node_type : str
            Specify the type of node the target is. Defaults to "protein", but can be any BEL node type.
        pmods : list
            If the target node is a protein, protein modifications can be added to restrict searching to protein nodes
            with both the given name AND the given pmods. Can use BEL shorthand pmods (e.g. "pho", "ace", "me0") or
            pmod values (e.g. "protein phosphorylation").
        neighbor_edge_type : str
            The class of edges to consider when looking for interactors. Defaults to "E" which encompasses all edges
            types, but more specific classes can be given such as only BEL edge types ("bel"), causal edges "causal",
            or "increases"/"directly_increases".
        print_sql : bool
            If True, will print each query as it is executed in OrientDB.
        """
        populate()  # Populate tables if empty

        self.node_name = node_name.lower()
        self.node_type = node_type
        self.pmods = pmods
        self.edge = neighbor_edge_type
        self.edge_filters = " AND " + " AND ".join(neighbor_edge_filters) if neighbor_edge_filters else ""
        self.print_sql = print_sql

        self.results = None

    def __len__(self):
        return len(self.results)

    def __query_graphstore(self, sql: str) -> Optional[pd.DataFrame]:
        if self.print_sql:
            print(sql)

        results = rest_query.sql(sql).table

        if isinstance(results, str):
            logger.warning("No results!")
            return None

        else:
            return results

    def __query_db(self, druggable: bool = False):
        """Check if query results are stored in cache."""
        # TODO remove if not used again soon
        target = self.node_name.upper()  # Target symbol upper case for humans
        if druggable:
            rels = EDGE_MAPPER['causal']

        elif self.edge in EDGE_MAPPER:
            rels = EDGE_MAPPER[self.edge]

        else:
            rels = [self.edge]

        table = Druggable if druggable else General

        query = session().query(table).filter_by(target_symbol=target, target_type=self.node_type)
        results = pd.read_sql(query.statement, query.session.bind)

        # Filter by edges and pmods
        filtered_df = results[results['relation_type'].isin(rels)]
        if self.pmods:
            filtered_df = filtered_df[filtered_df['pmod_type'].isin(self.pmods)]

        return filtered_df if not filtered_df.empty else None

    def find_interactors(self) -> pd.DataFrame:
        """Return interactors of the target."""
        sql = INTERACTOR_QUERY

        cols = ["target_species", "pmid", "pmc", "interactor_type", "interactor_name", "interactor_bel",
                "relation_type", "target_bel", "target_type", "target_symbol", "pmod_type"]

        if self.node_type != 'protein' or not self.pmods:
            sql = sql.replace('MATCH {{class:pmod, as:pmod{}}}<-has__pmod-', 'MATCH')
            formatted_sql = sql.format(self.node_type, self.node_name, self.edge, self.edge_filters)

        else:
            if 'all' in self.pmods:
                pmod_condition = "type != '' or name != ''"

            else:
                pmod_condition = f"type in {self.pmods}"

            pmod_string = f", WHERE:({pmod_condition})"

            if 'pho' in self.pmods or 'all' in self.pmods:
                pmod_string = pmod_string.replace(")", " OR name like '%phosphorylat%')")

            formatted_sql = sql.format(pmod_string, self.node_type, self.node_name, self.edge, self.edge_filters)

        df_results = self.__query_graphstore(formatted_sql)

        if df_results is not None:
            self.results = df_results[cols]
            self.results['target_species'] = self.results['target_species'].fillna(0).astype(int)

        return self.results

    def druggable_interactors(self) -> pd.DataFrame:
        """Return all druggable interactors of the target.

        Requires specialized queries and therefore is separate from `find_interactors`.
        """
        cols = ['drug', 'capsule_interactor_type', 'capsule_interactor_bel', 'interactor_bel', 'interactor_type',
                'interactor_name', 'relation_type', 'target_bel', 'target_symbol', 'target_type',
                'pmid', 'pmc', 'rel_pub_year', 'rel_rid', 'drug_rel_rid', 'drug_rel_actions',
                'drugbank_id', 'chembl_id', 'pubchem_id', 'pmod_type']

        if self.node_type != 'protein' or not self.pmods:
            pure_query = PURE_DRUGGABLE_QUERY.replace('MATCH {{class:pmod, as:pmod{}}}<-has__pmod-', 'MATCH')
            capsule_query_1 = CAPSULE_DRUGGABLE_MODIFIED.replace(
                'MATCH {{class:pmod, as:pmod{}}}<-has__pmod-', 'MATCH'
            )
            capsule_query_2 = CAPSULE_DRUGGABLE_COMPLEX.replace(
                'MATCH {{class:pmod, as:pmod{}}}<-has__pmod-', 'MATCH'
            )
            formatted_pure_sql = pure_query.format(self.node_type, self.node_name, self.edge_filters)
            formatted_capsule_sql_1 = capsule_query_1.format(self.node_type, self.node_name, self.edge_filters)
            formatted_capsule_sql_2 = capsule_query_2.format(self.node_type, self.node_name, self.edge_filters)

        else:
            if 'all' in self.pmods:
                pmod_condition = "type != '' or name != ''"
            else:
                pmod_condition = f"type in {self.pmods}"

            pmod_string = f", WHERE:({pmod_condition})"

            if 'pho' in self.pmods or 'all' in self.pmods:
                pmod_string = pmod_string.replace(")", " OR name like '%phosphorylat%')")

            # Drugs only for humans so only check one
            formatted_pure_sql = PURE_DRUGGABLE_QUERY.format(
                pmod_string, self.node_type, self.node_name, self.edge_filters
            )
            formatted_capsule_sql_1 = CAPSULE_DRUGGABLE_MODIFIED.format(
                pmod_string, self.node_type, self.node_name, self.edge_filters
            )
            formatted_capsule_sql_2 = CAPSULE_DRUGGABLE_COMPLEX.format(
                pmod_string, self.node_type, self.node_name, self.edge_filters
            )

        logger.info("Querying database...")

        pure_results = self.__query_graphstore(sql=formatted_pure_sql)
        capsule_results_1 = self.__query_graphstore(sql=formatted_capsule_sql_1)
        capsule_results_2 = self.__query_graphstore(sql=formatted_capsule_sql_2)

        results_check = [x is not None for x in (pure_results, capsule_results_1, capsule_results_2)]

        if any(results_check):  # Need only 1 not to be None
            df_concat = pd.concat(
                [pure_results, capsule_results_1, capsule_results_2], axis=0
            ).reindex(columns=cols)
            self.results = df_concat[cols]
            self.results["drug_rel_actions"] = self.results["drug_rel_actions"].str.join("|")
            self.results = self.results.drop_duplicates()

        return self.results

    def drug_and_interactors(self):
        """Return a list of interactors and the drugs that affect them."""
        if self.results is not None and 'drug' in self.results.columns:
            return self.results[['drug', 'interactor_name']]

    def unique_interactors(self):
        """Return a list of unique interactors found in the results dataframe."""
        if self.results is not None:
            return tuple(self.results['interactor_name'].unique())

    def unique_drugs(self):
        """Return a list of unique drugs found in the results dataframe."""
        if self.results is not None:
            return tuple(self.results['drug'].unique())


def get_interactor_list(results_df: pd.DataFrame):
    """Collect interactors from InteractorFinder results."""
    interactors = set()
    for gene_list in results_df.interactor_involved_genes:
        if gene_list is not None:
            for gene in gene_list:
                interactors.add(gene)
    for other_list in results_df.interactor_involved_other:
        if other_list is not None:
            for other in other_list:
                interactors.add(other)
    for name in results_df.interactor_name:
        interactors.add(name)

    return interactors
