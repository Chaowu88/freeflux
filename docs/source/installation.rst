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
  
Solver installation
-------------------
 
FreeFlux requires the numerical optimization framework `OpenOpt <https://openopt.org/>`_ for nonlinear regression. It can be installed by:
 
.. code-block:: python
  
  pip install openopt
  pip install FuncDesigner
  
FreeFlux utilizes the modeling language Pyomo to formulate linear optimization problem. By default, solver is not installed together with Pyomo, and thus should be installed independently. For example, to install glpk

.. code-block:: python
  
  conda install -c conda-forge glpk  
  
Dependencies and Compatibility 
------------------------------

The OpenOpt framework works well with Python 3.7, but may have problem with 3.8 and above. The function *clock()* in Python built-in module `time` was removed since Python 3.8, so manual correction of the installed openopt package is needed for compatible use. Specifically, ``clock`` in either import statement or function calls should be replaced with ``perf_counter`` in scripts *ooIter.py*, *runProbSolver.py* and *result.py*. Alternatively, one can use the `corrected ones <https://github.com/Chaowu88/freeflux/tree/main/openopt_patch>`__ to overwrite those with the same name in the installation path.
  
.. Note::
  If exception "ModuleNotFoundError: No module named 'numpy'" is raised during the installation of openopt or FuncDesigner, please install numpy first.

.. code-block:: python

  pip install "numpy>=1.20,<1.23"
  

