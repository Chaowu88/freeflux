Installation
============

Using PIP
---------

FreeFlux was tested in Python 3.7 and 3.8, it can be installed with *pip* from PyPI:

.. code-block:: python

  python -m pip install --upgrade pip
  pip install freeflux

or from source (install `git <https://git-scm.com/>`__ first):

.. code-block:: python

  git clone https://github.com/Chaowu88/freeflux.git /path/to/freeflux
  pip install /path/to/freeflux

.. Note::
  Installation within a `virtual environment <https://docs.python.org/3.8/tutorial/venv.html>`_ is recommended.
  
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

The OpenOpt framework works well with Python 3.7, but may have problem with 3.8 and above. The function *clock()* in Python built-in module `time` was removed since Python 3.8, so manual correction of the installed openopt package is needed for compatible use. Specifically, ``clock`` in either import statement or function calls should be replaced with ``perf_counter`` in scripts *ooIter.py*, *runProbSolver.py* and *result.py*.
  
.. Note::
  If exception "ModuleNotFoundError: No module named 'numpy'" is raised during the installation of openopt or FuncDesigner, please install numpy first.

.. code-block:: python

  pip install numpy~=1.20.2
  

