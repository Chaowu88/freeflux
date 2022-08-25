Installation
============

Using PIP
---------
FreeFlux was developed and tested in Python 3.8, it can be installed with *pip* from PyPI:
::
  python -m pip install --upgrade pip
  pip install freeflux
or from source:
::
  git clone https://github.com/Chaowu88/freeflux.git /path/to/freeflux
  pip install /path/to/freeflux
.. Note::
  Installation within a `virtual environment <https://docs.python.org/3.8/tutorial/venv.html>`_ is recommendated.
  
Additional Dependencies
-----------------------
FreeFlux also requires the numerical optimization framework `OpenOpt <https://openopt.org/>`_ for nonlinear regression. This framework can be installed using the following command:
::
  pip install openopt == 0.5629
  pip install FuncDesigner == 0.5629
.. Warning::
  The function *clock()* in Python built-in module time was removed since version 3.8, so manual correction of the installed openopt package is needed for compatible use. Specifically, *clock* in either import statement or function calls should be replaced with *perf_counter* in scripts *ooIter.py*, *runProbSolver.py* and *result.py*.
  
