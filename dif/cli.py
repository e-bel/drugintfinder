"""Console script for dif."""
import sys
import click

from dif.finder import InteractorFinder


@click.group(help=f"Druggable Interactor Finder Framework Command Line Utilities on {sys.executable}")
@click.version_option()
def main():
    """Console script for dif."""
    pass


@main.command()
@click.argument('symbol')
@click.option('-n', '--node', default='protein', help="Target node type. Defaults to 'protein'.")
@click.option('-e', '--edge', default='causal', help="Interactor/target relationship type. Defaults to 'causal'.")
@click.option('-p', '--pmods', default=[], help="Comma separated list of acceptable target protein modifications.")
@click.option('-d', '--druggable', is_flag=True, default=False, help="Flag to enable filtering of only druggable ints.")
@click.option('-s', '--sql', is_flag=True, default=False, help="Flag to print query.")
def find(symbol: str, node: str, edge: str, pmods: str, druggable: bool, sql: bool):
    """Identifies interactors of given target and criteria.

    Parameters
    ----------
    symbol: str
        Gene symbol of the target node.
    node: str
        Type of target node to consider (e.g. 'protein', 'rna', 'gene', etc.)
    edge: str
        Edge type between interactor and target nodes (e.g. 'increases', 'causal', 'correlative', 'E' for all, etc.)
    pmods: str
        Comma separated list of 3-letter target protein modifications. This will filter target node results for those
        with specified pmods.
    druggable: bool
        If True, will only include interactors that are targeted by a drug in the graph.
    sql: bool
        Flag to print query.

    Returns
    -------

    """
    if isinstance(pmods, str):
        pmods = pmods.split(",")

    finder = InteractorFinder(symbol=symbol, pmods=pmods, edge=edge)
    if not druggable:
        results = finder.find_interactors(target_type=node, print_sql=sql)

    else:
        results = finder.druggable_interactors(target_type=node, print_sql=sql)

    click.echo(results)

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
