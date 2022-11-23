FreeFlux
========

FreeFlux is a Python package for :sup:`13`\ C metabolic flux analysis of biological systems at isotopic steady state or transient state which thus can be used for both heterotrophic and autotrophic organisms. Functionally, FreeFlux is capable of:

- metabolic flux estimation
- simulation of labeling patterns of metabolite (fragments)
- constraint-based optimizations, such as flux balance analysis and flux variability analysis

The aim is to benefit the accessibility of :sup:`13`\ C fluxomics technique for researchers in the community of metabolic phenotyping and engineering.

A documentation can be found `here <https://freeflux.readthedocs.io/en/latest/index.html>`_. The documentation illustrates the fundamental functions of FreeFlux with a `toy model <https://github.com/Chaowu88/freeflux/tree/main/models/toy>`_. Two practical models of `E. coli <https://github.com/Chaowu88/freeflux/tree/main/models/ecoli>`_ and `Synechocystis <https://github.com/Chaowu88/freeflux/tree/main/models/synechocystis>`_ are also provided with `tutorials <https://github.com/Chaowu88/freeflux/tree/main/tutorials>`_.

Installation
============

FreeFlux was tested in Python 3.7 and 3.8. It can be installed using *pip* from PyPI:

.. code-block:: python

  python -m pip install --upgrade pip
  pip install freeflux

or from source:

.. code-block:: python

  git clone https://github.com/Chaowu88/freeflux.git /path/to/freeflux
  pip install /path/to/freeflux
  
Installation in an `virtual environment <https://docs.python.org/3.8/tutorial/venv.html>`_ is recommendated.

.. Note::
  FreeFlux requires the numerical optimization framework `OpenOpt <https://openopt.org/>`_ for nonlinear regression. It works well in Python 3.7 but may has compatibility issues in Python 3.8 and above. Please see `here <https://freeflux.readthedocs.io/en/latest/installation.html#dependency-compatibility>`_ for solutions.

Example Usage
=============

A typical use of FreeFlux starts with building a model by reading metabolic reactions with atom transitions. The model can call a handler of fitter, simulator or optimizor to perform flux estimation, labeling pattern simulation and constraint-based flux analysis, respectively. Different methods are provided for these handlers for data input and computation. 
Here is an example script of flux estimation at steady state using the `toy model <https://github.com/Chaowu88/freeflux/tree/main/models/toy>`_.

.. code-block:: python
   
   MODEL_FILE = 'path/to/reactions.tsv'
   MEASURED_MDVS = 'path/to/measured_MDVs.tsv'
   MEASURED_FLUXES = 'path/to/measured_fluxes.tsv'
   
   from freeflux import Model
   
   model = Model('demo')
   model.read_from_file(MODEL_FILE)
   
   with model.fitter('ss') as fit:
       fit.set_labeling_strategy('AcCoA', 
                                 labeling_pattern = ['01', '11'], 
                                 percentage = [0.25, 0.25], 
                                 purity = [1, 1])
       fit.set_flux_bounds('all', bounds = [-100, 100])
       fit.set_measured_MDV('Glu_12345', 
                            mean = [0.328,0.276,0.274,0.088,0.03,0.004], 
                            sd = [0.01,0.01,0.01,0.01,0.01,0.01])
       fit.set_measured_flux('v1', mean = 10, sd = 1)
       fit.prepare()
       res = fit.solve()

For more information, please see the `documentation <https://freeflux.readthedocs.io/en/latest/index.html>`.
