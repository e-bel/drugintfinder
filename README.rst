Druggable Interactor Finder |build| |docs| |zenodo|
==============================================================

A pipeline for identifying the most notable druggable interactors of a target within a e(BE:L) generated
Knowledge Graph. This feature will be integrated into the e(BE:L) package in the future.


Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------

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

.. |docs| image:: https://readthedocs.org/projects/druggable-interactor-finder/badge/?version=latest
        :target: https://druggable-interactor-finder.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. |pypi_license| image:: https://img.shields.io/pypi/l/drugintfinder.svg
    :target: https://pypi.python.org/pypi/drugintfinder
    :alt: MIT License

.. |pypi_version| image:: https://img.shields.io/pypi/v/drugintfinder.svg
    :alt: Current version on PyPI

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/drugintfinder.svg
    :alt: Stable Supported Python Versions

.. |coverage| image:: https://codecov.io/gh/e-bel/drugintfinder/coverage.svg?branch=master
    :target: https://codecov.io/gh/e-bel/drugintfinder?branch=master
    :alt: Coverage Status

.. |build| image:: https://travis-ci.com/e-bel/drugintfinder.svg?branch=master
    :target: https://travis-ci.com/e-bel/drugintfinder
    :alt: Build Status
    
.. |zenodo| image:: https://zenodo.org/badge/372590270.svg
   :target: https://zenodo.org/badge/latestdoi/372590270
