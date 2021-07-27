Druggable Interactor Finder's Documentation
==============================================================

A pipeline for identifying the most notable druggable interactors of a target within a e(BE:L) generated
Knowledge Graph.


Installation
--------------

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


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   cli
   installation
   usage
