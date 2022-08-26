Quick Start
===========

To begin with, we show here how to solve the metabolic fluxes at steady state with a toy model which includes the typical reactions in the TCA cycle with acetyl-CoA and aspartate as the initial substrates.

Loading a Toy Model
-------------------

A model can be defined in a tab-separated values (.tsv) or Excel spreadsheet (.xlsx) file with the format as below.

.. list-table:: reactions.tsv
   :widths: 25 50 50 15
   :header-rows: 1

   * - #reaction_ID
     - reactant_IDs(atom)
     - product_IDs(atom)
     - reversibility
   * - v1
     - OAA(abcd)+AcCoA(ef)
     - Cit(dcbfea)
     - 0
   * - v2
     - Cit(abcdef)
     - AKG(abcde)+CO2(f)
     - 0
   * - v3
     - AKG(abcde)
     - Glu(abcde)
     - 0
   * - v4
     - AKG(abcde)
     - Suc(bcde)+CO2(a)
     - 0
   * - v5
     - Suc(abcd,dcba)
     - Fum(abcd,dcba)
     - 0
   * - v6
     - Fum(abcd,dcba)
     - OAA(abcd)
     - 1
   * - v7
     - Asp(abcd)
     - OAA(abcd)
     - 0
     
.. Note::
  1. Letters in parenthesis indicates the C atom mapping in reactions.
  2. ``#`` is required in the header.
  
The model can then be loaded by:

.. code-block:: python
   
   MODEL_FILE = '../models/toy/reactions.tsv'
   
   from freeflux import Model
   
   model = Model('demo')
   model.read_from_file(MODEL_FILE)
   
Specifying the Labeling Strategy
--------------------------------

The small metabolic network uptakes 25% (mol%) C2 labeled acetyl-CoA, [2-\ :sup:`13`\C] AcCoA and 25% fully labeled acetyl-CoA, [U-\ :sup:`13`\C\ :sub:`2`\] AcCoA assuming 100% purity of the tracers, respectively. This labeling strategy can be set with the following line:

.. code-block:: python
  
  fit = model.fitter('ss')
  fit.set_labeling_strategy('AcCoA', ['01', '11'], [0.25, 0.25], [1, 1])

Adding Bounds for Fluxes
------------------------

The search range for all fluxes or a specific flux can be set using ``set_flux_bounds``. Here we constrained fluxes in the range of -100 and 100.

.. code-block:: python

  fit.set_flux_bounds('all', bounds = [-100, 100])
  
Reading the Measuremnets
------------------------

