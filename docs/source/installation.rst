Installation
============

Using PIP
---------

FreeFlux is compatible with Python versions 3.7 through 3.10 and can be easily installed using *pip* from PyPI. To get started, first upgrade pip using the following command:

.. code-block:: python

  python -m pip install --upgrade pip

Next, install FreeFlux with the following command:

.. code-block:: python

  pip install freeflux  

Alternatively, you can install FreeFlux from the source code by cloning the GitHub repository using the following command (assuming you have `git <https://git-scm.com/>`__ installed):

.. code-block:: python

  git clone https://github.com/Chaowu88/freeflux.git /path/to/freeflux

Then, install FreeFlux using *pip*:

.. code-block:: python

  pip install /path/to/freeflux
  
.. Note::
  Note that it's recommended to install FreeFlux within a `virtual environment <https://docs.python.org/3.8/tutorial/venv.html>`_ to avoid conflicts with other Python packages.

If FreeFlux is installed on a Linux platform with Python>=3.8, using `JAX <https://github.com/google/jax>`__ can significantly speed up the ``prepare`` step before flux estimation. JAX can be installed using the following command:

.. code-block:: python

  pip install "jax[cpu]>=0.4,<=0.4.18"

Solver Installation
-------------------
 
FreeFlux requires the numerical optimization framework `OpenOpt <https://openopt.org/>`_ for nonlinear regression, which can be installed with:
 
.. code-block:: python
  
  pip install openopt
  pip install FuncDesigner
  
In addition, FreeFlux uses the modeling language `Pyomo <http://www.pyomo.org/>`__ to formulate linear optimization problems. If the solvers are not installed together with Pyomo, you will need to install them independently. For example, to install the glpk solver:

.. code-block:: python
  
  conda install -c conda-forge glpk  
  
Dependencies and Compatibility 
------------------------------

The OpenOpt framework works well with Python 3.7, but may have compatibility issues with 3.8 and above due to the removal of the *clock()* function in Python's built-in module `time`. To use FreeFlux with Python 3.8 and above, manual correction of the installed OpenOpt package is needed. Specifically, ``clock`` in either import statement or function calls should be replaced with ``perf_counter`` in scripts ooIter.py, runProbSolver.py and result.py. Alternatively, one can use the `corrected ones <https://github.com/Chaowu88/freeflux/tree/main/openopt_patch>`__ to overwrite those with the same name in the installation path.
  
.. Note::
  Note: If you encounter the "ModuleNotFoundError: No module named 'numpy'" exception during the installation of OpenOpt or FuncDesigner, please install numpy first using the following command:

.. code-block:: python

  pip install "numpy>=1.20,<1.23"
  

