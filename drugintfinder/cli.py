"""Console script for drugintfinder."""
import sys

import click

from drugintfinder.ranker import Ranker
from drugintfinder.utils import export_table
from drugintfinder.finder import InteractorFinder


@click.group(help=f"Druggable Interactor Finder Framework Command Line Utilities on {sys.executable}")
@click.version_option()
def main():
    """Console script for drugintfinder."""
    pass


@main.command()
@click.argument('name')
@click.option('-n', '--node-type', default='protein', help="Target node type. Defaults to 'protein'.")
@click.option('-e', '--edge-type', default='causal', help="Interactor/target relationship type. Defaults to 'causal'.")
@click.option('-m', '--pmods', default=[], help="Comma separated list of acceptable target protein modifications.")
@click.option('-d', '--druggable', is_flag=True, default=False, help="Flag to enable filtering of druggable ints.")
@click.option('-s', '--sql', is_flag=True, default=False, help="Flag to print query.")
@click.option('-o', '--output', default=None, help="Results output path - defaults to Excel.")
@click.option('-v', '--verbose', is_flag=True, default=False, help="Flag to print results to STDOUT.")
def find(
        name: str, node_type: str, edge_type: str, pmods: str, druggable: bool, sql: bool, output: str, verbose: bool
):
    """Identify interactors of given target and criteria.

    Parameters
    ----------
    name: str
        Gene symbol or namespace value of the target node.
    node_type: str
        Type of target node to consider (e.g. 'protein', 'rna', 'gene', etc.)
    edge_type: str
        Edge type between interactor and target nodes (e.g. 'increases', 'causal', 'correlative', 'E' for all, etc.)
    pmods: str
        Comma separated list of 3-letter target protein modifications. This will filter target node results for those
        with specified pmods.
    druggable: bool
        If True, will only include interactors that are targeted by a drug in the graph.
    sql: bool
        Flag to print query.
    output: str
        Path to write results to. Currently only supports Excel.
    verbose: bool
        Flag to print results table to STDOUT.
    """
    if isinstance(pmods, str):
        pmods = pmods.split(",")

    finder = InteractorFinder(
        node_name=name, pmods=pmods, neighbor_edge_type=edge_type, node_type=node_type, print_sql=sql
    )
    if not druggable:
        finder.find_interactors()

    else:
        finder.druggable_interactors()

    if output:
        export_table(finder.results, output)

    if verbose:
        click.echo(finder.results)


@main.command()
@click.argument('name')
@click.option('-n', '--node-type', default='protein', help="Target node type. Defaults to 'protein'.")
@click.option('-m', '--pmods', default=[], help="Comma separated list of acceptable target protein modifications.")
@click.option('-r', '--reward', default=1, help="Points awarded for passing inspection criteria.")
@click.option('-p', '--penalty', default=-1, help="Points penalized for failing inspection criteria.")
@click.option('-o', '--output', default=None, help="Results output path - defaults to Excel.")
@click.option('-t', '--pivot', is_flag=True, default=False, help="Flag to enable target focused summary results.")
@click.option('-v', '--verbose', is_flag=True, default=False, help="Flag to print results to STDOUT.")
def rank(name: str, node_type: str, pmods: str, reward: str, penalty: str, output: str, verbose: bool, pivot: bool):
    """Ranks the drug/interactor combos with metadata and returns a summary table.

    Parameters
    ----------
    name: str
        Gene symbol or namespace value of the target node.
    node_type: str
        Type of target node to consider (e.g. 'protein', 'rna', 'gene', etc.)
    pmods: str
        Comma separated list of 3-letter target protein modifications. This will filter target node results for those
        with specified pmods.
    reward: str
        Number of points to reward entries that follow desired outcome.
    penalty: str
        Number of points to penalize entries that follow desired outcome.
    output: str
        Path to write results to. Currently only supports Excel.
    pivot: bool
        Flag to enable target focused summary results.
    verbose: bool
        Flag to print results table to STDOUT.
    """
    if isinstance(pmods, str):
        pmods = pmods.split(",")

    ranker = Ranker(name=name, node_type=node_type, pmods=pmods, reward=int(reward), penalty=int(penalty))
    ranker.rank()
    summary = ranker.summarize(pivot=pivot)

    if output:
        export_table(summary, output)

    if verbose:
        click.echo(summary)


if __name__ == "__main__":
    sys.exit(main())
