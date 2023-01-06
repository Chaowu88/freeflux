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
   
   MODEL_FILE = 'path/to/reactions.tsv'
   MEASURED_MDVS = 'path/to/measured_MDVs.tsv'
   MEASURED_FLUXES = 'path/to/measured_fluxes.tsv'
   
   from freeflux import Model
   
   model = Model('demo')
   model.read_from_file(MODEL_FILE)
   
Specifying the Labeling Strategy
--------------------------------

The small metabolic network uptakes 25% (mol%) C2 labeled acetyl-CoA, [2-\ :sup:`13`\C] AcCoA and 25% fully labeled acetyl-CoA, [U-\ :sup:`13`\C\ :sub:`2`\] AcCoA assuming 100% purity of these two tracers, respectively. The labeling strategy can be set with the following line:

.. code-block:: python
  
  fit = model.fitter('ss')
  fit.set_labeling_strategy('AcCoA', ['01', '11'], [0.25, 0.25], [1, 1])

.. Note::
   1. If the sum of percantage of the specified isotopomer tracers is less than 1, the remaining 1-sum will be considered as the unlabeled form (atoms in natural abundance).
   2. Call this method for each substrate if multiple labeled substrates are used.
   
Adding Bounds for Fluxes
------------------------

The search range for all fluxes or a specific flux can be set using ``set_flux_bounds``. Here we constrain the fluxes to the range from -100 to 100.

.. code-block:: python

  fit.set_flux_bounds('all', bounds = [-100, 100])
  
Reading the Measuremnets
------------------------

The labeling pattern or mass isotopomer distribution vector (MDV) of measurable metabolites or metabolite fragments should be provided for the fitting as well as measured exchange reactions, e.g., substrate uptake, product secretion and specific growth rate. They can be provided using the lines:

.. code-block:: python
   
   fit.set_measured_MDV('Glu_12345', 
                        mean = [0.328,0.276,0.274,0.088,0.03,0.004], 
                        sd = [0.01,0.01,0.01,0.01,0.01,0.01])
   fit.set_measured_flux('v1', mean = 10, sd = 1)

.. Note::
   For input of a set of measured MDVs and fluxes, it is more convenient to read from a .tsv or xlsx. file by methods ``set_measured_MDVs_from_file`` and ``set_measured_fluxes_from_file``.
   
Solve the Fluxes
----------------

Now we can solve the flux distribution in the toy model by:

.. code-block:: python
   
   fit.prepare()
   res = fit.solve()
   
The ``solve`` method returns a FitResults object. The estimated net and total (including both forward and backward fluxes in reversible reactions) fluxes can be accessed by the attributes ``opt_net_fluxes`` and ``opt_total_fluxes``.

With Statement
--------------

The returned Fitter object by calling ``fitter`` method is actually a context manager, and thus the above flux estimation can also be done using the with statement:

.. code-block:: python
   
   with model.fitter('ss') as fit:
       fit.set_labeling_strategy('AcCoA', ['01', '11'], [0.25, 0.25], [1, 1])
       fit.set_flux_bounds('all', bounds = [-100, 100])
       fit.set_measured_MDV('Glu_12345', 
                            mean = [0.328,0.276,0.274,0.088,0.03,0.004], 
                            sd = [0.01,0.01,0.01,0.01,0.01,0.01])
       fit.set_measured_flux('v1', mean = 10, sd = 1)
       fit.prepare()
       res = fit.solve()
       
       
