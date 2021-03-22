"""Methods for executing pipeline."""
import logging
import pandas as pd

from ebel_rest import query as rest_query

from dif.constants import INTERACTOR_QUERY, PURE_DRUGGABLE_QUERY, CAPSULE_DRUGGABLE_QUERY

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class InteractorFinder:

    def __init__(self, name_list: str):
        self.names = name_list
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

    def find_interactors(self,
                         target_type: str = 'protein',
                         edge_class: str = 'E',
                         pmods: list = None,
                         print_sql: bool = False) -> pd.DataFrame:
        """Returns causal interactors of the target

        Parameters
        ----------
        pmods : list
            A python list of target protein modifications for filtering results. Pmods must be defined using the
            e(BE:L) pmod labels.

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
        """Returns all interactors of the target

        Parameters
        ----------
        pmods : list
            A python list of target protein modifications for filtering results. Pmods must be defined using the
            e(BE:L) pmod labels.

        Returns
        -------
        pd.DataFrame
        """
        pure_query = PURE_DRUGGABLE_QUERY
        capsule_query = CAPSULE_DRUGGABLE_QUERY

        cols = ['drug', 'capsule_interactor_type', 'capsule_interactor_bel', 'interactor_bel',
                'interactor_type', 'interactor_name', 'relation_type', 'target_bel', 'pmid', 'pmc',
                'rel_rid', 'drug_rel_rid']

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
        """Exports results dataframe to path"""

        if self.results is None:
            logger.warning("No results found! Failed to export.")

        else:
            self.results.to_excel(file_path, index=False)
            logger.info(f"Results written to {file_path}")

    def drug_and_interactors(self):
        """Returns a list of interactors and the drugs that affect them"""
        if self.results is not None and 'drug' in self.results.columns:
            return self.results[['drug', 'interactor_name']]

    def unique_drugs(self):
        """Returns a list of unique drugs found in the results dataframe"""
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
