"""Methods for executing pipeline."""
import logging
import pandas as pd

from ebel_rest import query as rest_query

from dif.defaults import session
from dif.models import General, Druggable
from dif.constants import INTERACTOR_QUERY, PURE_DRUGGABLE_QUERY, CAPSULE_DRUGGABLE_QUERY, EDGE_MAPPER

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class InteractorFinder:

    def __init__(self, symbol: str):
        self.names = list({symbol, symbol.upper(), symbol.lower(), symbol.capitalize()})
        self.results = None

    def __len__(self):
        return len(self.results)

    @staticmethod
    def __query(sql: str, print_sql: bool = False) -> pd.DataFrame:
        if print_sql:
            print(sql)
        results = rest_query.sql(sql)
        if results:
            return results.table
        else:
            logger.warning("No results!")

    def __check_db(self, node_type: str, edge_type: str, pmods: list, druggable: bool = False):
        """Checks if query results are stored in cache."""
        target = self.names[0].upper()  # Target symbol upper case for humans
        if edge_type in EDGE_MAPPER:
            rels = EDGE_MAPPER[edge_type]

        else:
            rels = [edge_type]

        table = Druggable if druggable else General

        query = session().query(table).filter_by(target_symbol=target, target_type=node_type).all()
        results = pd.read_sql(query.statement, session().query.bind)

        # Filter by edges and pmods
        filtered_df = results[results['relation_type'].isin(rels)]
        if pmods:
            filtered_df = filtered_df[filtered_df['pmod_type'].isin(pmods)]

        return filtered_df if not filtered_df.empty else None


    def find_interactors(self,
                         target_type: str = 'protein',
                         edge_class: str = 'E',
                         pmods: list = None,
                         print_sql: bool = False) -> pd.DataFrame:
        """Returns interactors of the target.

        Parameters
        ----------
        target_type: str
            Node type to query.
        edge_class: str
            Which edge type to restrict the search to.
        pmods : list
            A python list of target protein modifications for filtering results. Pmods must be defined using the
            e(BE:L) pmod labels.
        print_sql: bool
            If true, will print the query.

        Returns
        -------
        pd.DataFrame
        """
        sql = INTERACTOR_QUERY

        cols = ["target_species", "pmid", "pmc", "interactor_type", "interactor_involved_genes",
                "interactor_involved_other", "interactor_name", "interactor_bel", "relation_type",
                "target_bel"]

        if target_type != 'protein':
            sql = sql.replace('MATCH {{class:pmod, as:pmod{}}}<-has__pmod-', 'MATCH')
            formatted_sql = sql.format(target_type, self.names, edge_class)

        elif pmods:
            cols.append("pmod_type")

            if 'all' in pmods:
                pmod_condition = "type != '' or name != ''"
            else:
                pmod_condition = f"type in {pmods}"

            pmod_string = f", WHERE:({pmod_condition})"

            if 'pho' in pmods or 'all' in pmods:
                pmod_string = pmod_string.replace(")", " OR name like '%phosphorylat%')")
            formatted_sql = sql.format(pmod_string, target_type, self.names, edge_class)

        else:
            formatted_sql = sql.format("", target_type, self.names, edge_class)

        df_results = self.__query(formatted_sql, print_sql=print_sql)

        self.results = df_results[cols]

        return self.results

    def druggable_interactors(self,
                              target_type: str = 'protein',
                              pmods: list = None,
                              print_sql: bool = False) -> pd.DataFrame:
        """Returns all druggable interactors of the target. Requires specialized queries and therefore is separate from
        `find_interactors`.

        Parameters
        ----------
        target_type: str
            Node type to query.
        pmods : list
            A python list of target protein modifications for filtering results. Pmods must be defined using the
            e(BE:L) pmod labels.
        print_sql: bool
            If true, will print the query.

        Returns
        -------
        pd.DataFrame
        """
        pure_query = PURE_DRUGGABLE_QUERY
        capsule_query = CAPSULE_DRUGGABLE_QUERY

        cols = ['drug', 'capsule_interactor_type', 'capsule_interactor_bel', 'interactor_bel', 'interactor_type',
                'interactor_name', 'relation_type', 'target_bel', 'pmid', 'pmc', 'rel_pub_year', 'rel_rid',
                'drug_rel_actions', 'drug_rel_rid', 'evidence', 'drugbank_id', 'drug_patents', 'drug_products',
                'chembl_id', 'pubchem_id']

        if target_type != 'protein':
            pure_query = pure_query.replace('MATCH {{class:pmod, as:pmod{}}}<-has__pmod-', 'MATCH')
            capsule_query = capsule_query.replace('MATCH {{class:pmod, as:pmod{}}}<-has__pmod-', 'MATCH')
            formatted_pure_sql = pure_query.format(target_type, self.names)
            formatted_capsule_sql = capsule_query.format(target_type, self.names)

        elif pmods:
            if 'all' in pmods:
                pmod_condition = "type != '' or name != ''"
            else:
                pmod_condition = f"type in {pmods}"

            pmod_string = f", WHERE:({pmod_condition})"

            if 'pho' in pmods or 'all' in pmods:
                pmod_string = pmod_string.replace(")", " OR name like '%phosphorylat%')")

            formatted_pure_sql = pure_query.format(pmod_string, target_type, self.names)
            formatted_capsule_sql = capsule_query.format(pmod_string, target_type, self.names)

        else:
            formatted_pure_sql = pure_query.format("", target_type, self.names)
            formatted_capsule_sql = capsule_query.format("", target_type, self.names)

        logger.info("Querying database...")

        pure_results = self.__query(sql=formatted_pure_sql, print_sql=print_sql)
        capsule_results = self.__query(sql=formatted_capsule_sql, print_sql=print_sql)

        df_concat = pd.concat([pure_results, capsule_results], axis=0)

        self.results = df_concat[cols]

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
