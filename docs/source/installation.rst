Installation
============

Using PIP
---------

FreeFlux was tested in Python 3.7 and 3.8, it can be installed with *pip* from PyPI:

.. code-block:: python

  python -m pip install --upgrade pip
  pip install freeflux

or from source:

.. code-block:: python

  git clone https://github.com/Chaowu88/freeflux.git /path/to/freeflux
  pip install /path/to/freeflux

.. Note::
  Installation within a `virtual environment <https://docs.python.org/3.8/tutorial/venv.html>`_ is recommendated.
  
Dependencies and Compatibility 
------------------------

FreeFlux requires the numerical optimization framework `OpenOpt <https://openopt.org/>`_ for nonlinear regression. This framework works well with Python 3.7, but will have problem with 3.8 and above. The function *clock()* in Python built-in module `time` was removed since Python 3.8, so manual correction of the installed openopt package is needed for compatible use. Specifically, ``clock`` in either import statement or function calls should be replaced with ``perf_counter`` in scripts *ooIter.py*, *runProbSolver.py* and *result.py*.
  
.. Note::
  If exception "ModuleNotFoundError: No module named 'numpy'" is raised during the installation of openopt or FuncDesigner, please install numpy first.

.. code-block:: python

  pip install numpy~=1.20.2
  
FreeFlux also utilizes the modeling language Pyomo to formulate linear optimization problem. By default, solvers are not installed together with Pyomo, and thus should be installed independently.

.. code-block:: python
  
  conda install -c conda-forge glpk
