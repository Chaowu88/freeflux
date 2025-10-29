========
FreeFlux
========

FreeFlux is a Python package designed for :sup:`13`\ C metabolic flux analysis of biological systems at isotopic steady state or transient state, making it suitable for both heterotrophic and autotrophic organisms. With FreeFlux, you can:

- Estimate metabolic fluxes 
- Simulate labeling patterns of metabolite (fragments)
- Conduct constraint-based optimizations such as flux balance analysis and flux variability analysis

Our goal is to increase the accessibility of :sup:`13`\ C fluxomics techniques to researchers in the metabolic phenotyping and engineering community.

To get started, check out our `documentation <https://freeflux.readthedocs.io/en/latest/index.html>`__, which provides an overview of FreeFlux's fundamental functions using a `toy model <https://github.com/Chaowu88/freeflux/tree/main/models/toy>`__. We've also included `tutorials <https://github.com/Chaowu88/freeflux/tree/main/tutorials>`__ for practical models of `E. coli <https://github.com/Chaowu88/freeflux/tree/main/models/ecoli>`__ and `Synechocystis <https://github.com/Chaowu88/freeflux/tree/main/models/synechocystis>`__.

FreeFlux now supports :sup:`15`\ N labeling data, which can be enabled by specifying nitrogen as the labeled atom when defining the labeling strategy.

Installation
============

FreeFlux was tested in Python 3.7, 3.8, 3.9 and 3.10. It can be installed using *pip* from PyPI:

.. code-block:: python

  python -m pip install --upgrade pip
  pip install freeflux

or from source (assuming you have `git <https://git-scm.com/>`__ installed):

.. code-block:: python

  git clone https://github.com/Chaowu88/freeflux.git /path/to/freeflux
  pip install /path/to/freeflux
  
Installation within an `virtual environment <https://docs.python.org/3.8/tutorial/venv.html>`__ is recommendated.

If FreeFlux is installed on a Linux platform with Python>=3.8, using `JAX <https://github.com/google/jax>`__ can significantly speed up the preparation step before flux estimation. JAX can be installed using the following command:

.. code-block:: python

  pip install "jax[cpu]>=0.4,<=0.4.18"

Solver installation
===================

FreeFlux requires the numerical optimization framework `OpenOpt <https://openopt.org/>`__ for nonlinear regression. It can be installed by the following commands:

.. code-block:: python

  pip install openopt
  pip install FuncDesigner

Note that OpenOpt is known to work well in Python 3.7, but may have compatibility issues in Python 3.8 and above. FreeFlux now automatically applies the necessary patches to ensure compatibility. If you prefer to apply the patch manually, please refer to this `link <https://freeflux.readthedocs.io/en/latest/installation.html#dependencies-and-compatibility>`__ for solutions.

FreeFlux uses the modeling language Pyomo to formulate linear optimization problem. By default, solvers are not installed together with Pyomo and should be installed independently. For example, to install the glpk solver, run the following command:

.. code-block:: python
  
  conda install -c conda-forge glpk

Example Usage
=============

A typical workflow with FreeFlux starts by building a model through reading metabolic reactions with atom transitions. Users can then call a handler, such as the fitter, simulator, or optimizer, to perform flux estimation, labeling pattern simulation, or constraint-based flux analysis, respectively. Various methods are provided for these handlers for data input and computation.

Below is an example script that performs flux estimation at steady state using the `toy model <https://github.com/Chaowu88/freeflux/tree/main/models/toy>`__:

.. code-block:: python
   
   MODEL_FILE = 'path/to/reactions.tsv'
   
   from freeflux import Model
   
   model = Model('demo')
   model.read_from_file(MODEL_FILE)
   
   with model.fitter('ss') as fit:
       fit.set_labeling_strategy(
           'AcCoA', 
           labeling_pattern = ['01', '11'], 
           percentage = [0.25, 0.25], 
           purity = [1, 1],
           label_atom = 'C'
       )
       fit.set_flux_bounds('all', bounds = [-100, 100])
       fit.set_measured_MDV(
           'Glu_12345', 
           mean = [0.328,0.276,0.274,0.088,0.03,0.004], 
           sd = [0.01,0.01,0.01,0.01,0.01,0.01]
       )
       fit.set_measured_flux('v1', mean = 10, sd = 1)
       fit.prepare()
       res = fit.solve()

For more information, please refer to the `documentation <https://freeflux.readthedocs.io/en/latest/index.html>`__.

Performance Tip
===============
NumPy's underlying BLAS/LAPACK implementation may automatically utilize multiple CPU cores to improve computational efficiency. It can introduce additional overhead and may cause high system CPU usage, especially when using the ``solve_with_confidence_intervals`` method for confidence interval estimation. In such cases, it is recommended to explicitly limit the number of threads. For example, if OpenBLAS is used:

.. code-block:: bash

  export OPENBLAS_NUM_THREADS=1

License
=======

FreeFlux is released under the GPL version 3 license, please see `here <https://github.com/Chaowu88/freeflux/blob/main/LICENSE>`__ for more details.

Citation
========

`Wu et al.. FreeFlux: A Python Package for Time-Efficient Isotopically Nonstationary Metabolic Flux Analysis. ACS Synthetic Biology. 2023. <https://pubs.acs.org/doi/full/10.1021/acssynbio.3c00265>`__

Feel free to provide feedback at chao.wu@nrel.gov or chaowu09@gmail.com.
