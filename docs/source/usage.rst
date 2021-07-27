===========================
Druggable Interactor Finder
===========================

A pipeline for identifying the most notable druggable interactors of a target within a e(BE:L) generated
Knowledge Graph.


Installation
------------------------------------------------------

.. code-block:: sh

    $ pip install drugintfinder


Usage
--------
Finding the direct causal interactors of the phosphorylated MAPT protein:

.. code-block:: sh

    $ dif find MAPT -n protein -e causal -m pho

Finding the direct causal interactors of the phosphorylated MAPT protein which are druggable:

.. code-block:: sh

    $ dif find MAPT -n protein -e causal -m pho -d

Creating a table of various parameters for each interactor by which to rank them:

.. code-block:: sh

    $ dif rank MAPT -m pho


Disclaimer
----------

The Druggable Interactor finder is a resource developed in an academic capacity funded by the
`MAVO project <https://www.scai.fraunhofer.de/en/business-research-areas/bioinformatics/projects.html>`_
and thus comes with no warranty or guarantee of maintenance or support.


.. |pypi| image:: https://img.shields.io/pypi/v/drugintfinder.svg
        :target: https://pypi.python.org/pypi/drugintfinder

.. |travis| image:: https://img.shields.io/travis/e-bel/drugintfinder.svg
        :target: https://travis-ci.org/cebel/drugintfinder

.. |docs| image:: https://readthedocs.org/projects/drugintfinder/badge/?version=latest
        :target: https://ebel-rest.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. |pypi_license| image:: https://img.shields.io/pypi/l/drugintfinder.svg
    :target: https://pypi.python.org/pypi/drugintfinder
    :alt: MIT License

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/drugintfinder.svg
    :alt: Stable Supported Python Versions

.. |coverage| image:: https://codecov.io/gh/e-bel/drugintfinder/coverage.svg?branch=master
    :target: https://codecov.io/gh/e-bel/drugintfinder?branch=master
    :alt: Coverage Status

.. |build| image:: https://travis-ci.com/e-bel/drugintfinder.svg?branch=master
    :target: https://travis-ci.com/e-bel/drugintfinder
    :alt: Build Status
