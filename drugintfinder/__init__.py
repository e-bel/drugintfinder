"""A pipeline for identifying the most notable druggable interactors of a target within a e(BE:L) Knowledge Graph.

This package is designed for identifying the direct interactors of a specified node within an e(BE:L) generated
Knowledge Graph. The starting node can be filtered by given protein modifications and the type of interactors
(e.g. gene, RNA, protein) can also be controlled. Additioanlly, this package also can identify those interactors
that are druggable and provide a table of useful statistics about each drug/target combination that can be used for
ranking these pairs.

To get started with the druggable interactor finder, you can install it using:
``pip install drugintfinder``
"""

from drugintfinder import defaults

__author__ = """Bruce Schultz"""
__email__ = 'bruce.schultz@scai.fraunhofer.de'
__version__ = '0.2.2'
