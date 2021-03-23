"""Methods for executing pipeline."""
import logging
import pandas as pd

from typing import Optional
from ebel_rest import query as rest_query

from dif.defaults import session, engine
from dif.models import General, Druggable
from dif.constants import INTERACTOR_QUERY, PURE_DRUGGABLE_QUERY, CAPSULE_DRUGGABLE_QUERY, EDGE_MAPPER

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class InteractorFinder:

    def __init__(self, symbol: str, pmods: list = None, edge: str = 'E'):
        self.names = list({symbol, symbol.upper(), symbol.lower(), symbol.capitalize()})
        self.pmods = pmods
        self.edge = edge
        self.results = None

    def __len__(self):
        return len(self.results)

    @staticmethod
    def __query_graphstore(sql: str, print_sql: bool = False) -> Optional[pd.DataFrame]:
        if print_sql:
            print(sql)
        results = rest_query.sql(sql).table

        if isinstance(results, str):
            logger.warning("No results!")
            return None

        else:
            return results

    def __query_db(self, node_type: str, druggable: bool = False):
        """Checks if query results are stored in cache."""
        target = self.names[0].upper()  # Target symbol upper case for humans
        if druggable:
            rels = EDGE_MAPPER['causal']

        elif self.edge in EDGE_MAPPER:
            rels = EDGE_MAPPER[self.edge]

        else:
            rels = [self.edge]

        table = Druggable if druggable else General

        query = session().query(table).filter_by(target_symbol=target, target_type=node_type)
        results = pd.read_sql(query.statement, query.session.bind)

        # Filter by edges and pmods
        filtered_df = results[results['relation_type'].isin(rels)]
        if self.pmods:
            filtered_df = filtered_df[filtered_df['pmod_type'].isin(self.pmods)]

        return filtered_df if not filtered_df.empty else None

    def find_interactors(self,
                         target_type: str = 'protein',
                         print_sql: bool = False) -> pd.DataFrame:
        """Returns interactors of the target.

        Parameters
        ----------
        target_type: str
            Node type to query.
        print_sql: bool
            If true, will print the query.

        Returns
        -------
        pd.DataFrame
        """
        cached_results = self.__query_db(node_type=target_type, druggable=False)
        if cached_results is not None and not cached_results.empty:
            self.results = cached_results

        else:
            sql = INTERACTOR_QUERY

            cols = ["target_species", "pmid", "pmc", "interactor_type", "interactor_name", "interactor_bel",
                    "relation_type", "target_bel", "target_type", "target_symbol", "pmod_type"]

            if target_type != 'protein' or not self.pmods:
                sql = sql.replace('MATCH {{class:pmod, as:pmod{}}}<-has__pmod-', 'MATCH')
                formatted_sql = sql.format(target_type, self.names, self.edge)

            else:
                if 'all' in self.pmods:
                    pmod_condition = "type != '' or name != ''"
                else:
                    pmod_condition = f"type in {self.pmods}"

                pmod_string = f", WHERE:({pmod_condition})"

                if 'pho' in self.pmods or 'all' in self.pmods:
                    pmod_string = pmod_string.replace(")", " OR name like '%phosphorylat%')")
                formatted_sql = sql.format(pmod_string, target_type, self.names, self.edge)

            df_results = self.__query_graphstore(formatted_sql, print_sql=print_sql)

            if df_results is not None and not df_results.empty:
                self.results = df_results[cols]
                self.results['target_species'] = self.results['target_species'].fillna(0).astype(int)
                self.results.index += 1
                self.results.index.rename('id', inplace=True)

                logger.info(f"Importing results for {self.names[0].upper()} into SQLite DB")
                self.results.to_sql('general', if_exists="replace", con=engine)

        return self.results

    def druggable_interactors(self,
                              target_type: str = 'protein',
                              print_sql: bool = False) -> pd.DataFrame:
        """Returns all druggable interactors of the target. Requires specialized queries and therefore is separate from
        `find_interactors`.

        Parameters
        ----------
        target_type: str
            Node type to query.
        print_sql: bool
            If true, will print the query.

        Returns
        -------
        pd.DataFrame
        """
        cached_results = self.__query_db(node_type=target_type, druggable=True)
        if cached_results:
            self.results = cached_results

        else:
            pure_query = PURE_DRUGGABLE_QUERY
            capsule_query = CAPSULE_DRUGGABLE_QUERY

            cols = ['drug', 'capsule_interactor_type', 'capsule_interactor_bel', 'interactor_bel', 'interactor_type',
                    'interactor_name', 'relation_type', 'target_bel', 'target_symbol', 'target_type',
                    'pmid', 'pmc', 'rel_pub_year', 'rel_rid', 'drug_rel_actions', 'drug_rel_rid', 'evidence',
                    'drugbank_id', 'drug_patents', 'drug_products', 'chembl_id', 'pubchem_id', 'pmod_type']

            if target_type != 'protein' or not self.pmods:
                pure_query = pure_query.replace('MATCH {{class:pmod, as:pmod{}}}<-has__pmod-', 'MATCH')
                capsule_query = capsule_query.replace('MATCH {{class:pmod, as:pmod{}}}<-has__pmod-', 'MATCH')
                formatted_pure_sql = pure_query.format(target_type, self.names)
                formatted_capsule_sql = capsule_query.format(target_type, self.names)

            else:
                if 'all' in self.pmods:
                    pmod_condition = "type != '' or name != ''"
                else:
                    pmod_condition = f"type in {self.pmods}"

                pmod_string = f", WHERE:({pmod_condition})"

                if 'pho' in self.pmods or 'all' in self.pmods:
                    pmod_string = pmod_string.replace(")", " OR name like '%phosphorylat%')")

                # Drugs only for humans so only check one
                formatted_pure_sql = pure_query.format(pmod_string, target_type, [self.names[0].upper()])
                formatted_capsule_sql = capsule_query.format(pmod_string, target_type, [self.names[0].upper()])

            logger.info("Querying database...")

            pure_results = self.__query_graphstore(sql=formatted_pure_sql, print_sql=print_sql)
            capsule_results = self.__query_graphstore(sql=formatted_capsule_sql, print_sql=print_sql)

            df_concat = pd.concat([pure_results, capsule_results], axis=0)

            if df_concat is not None and not df_concat.empty:
                self.results = df_concat[cols]
                self.results.index += 1
                self.results.index.rename('id', inplace=True)

                logger.info(f"Importing results for {self.names[0].upper()} into SQLite DB")
                self.results.to_sql('druggable', if_exists="replace", con=engine)

        return self.results

    def export(self, file_path: str):
        """Exports results dataframe to path."""

        if self.results is None:
            logger.warning("No results found! Failed to export.")

        else:
            self.results.to_excel(file_path, index=False)
            logger.info(f"Results written to {file_path}")

    def drug_and_interactors(self):
        """Returns a list of interactors and the drugs that affect them."""
        if self.results is not None and 'drug' in self.results.columns:
            return self.results[['drug', 'interactor_name']]

    def unique_drugs(self):
        """Returns a list of unique drugs found in the results dataframe."""
        if self.results is not None:
            return pd.DataFrame(self.results['drug'].unique(), columns=['drug'])


def get_interactor_list(results_df: pd.DataFrame):
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
