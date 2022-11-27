# Configuration file for the Sphinx documentation builder.

# ------------------------ Project information ------------------------
project = 'freeflux'
copyright = '2022, Chao Wu'
author = 'Chao Wu'

import sys
from os.path import dirname, join

SRC_PATH = join(dirname(dirname(dirname(__file__))), 'src')
sys.path.insert(0, SRC_PATH)

from freeflux import __version__ as version

release = version
# version = '0.3.0'
# release = '0.3.0'


# ------------------------ General configuration ------------------------
extensions = ['sphinx.ext.duration',
              'sphinx.ext.doctest',
              'sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.intersphinx',
              'nbsphinx',
              'autoapi.extension']

autoapi_type = 'python'
autoapi_dirs = [join(SRC_PATH, 'freeflux')]

master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']

pygments_style = 'sphinx'






# ------------------------ Options for HTML output ------------------------
import sphinx_rtd_theme

html_theme = 'sphinx_rtd_theme'


# ------------------------ Options for Texinfo output ------------------------
# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [('index',
                      project,
                      project+' Documentation',
                      author,
                      'feeflux',
                      'A package for 13C metabolic flux analysis',
                      'Miscellaneous')]

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'python': ('https://docs.python.org/3/', None),
                       'numpy': ('https://numpy.org/doc/stable/', None),
                       'scipy': ('https://docs.scipy.org/doc/scipy/', None),
                       'openopt': ('https://openopt.org/Doc', None),
                       'pyomo': ('https://pyomo.readthedocs.io/en/stable/', None),
                       'sympy': ('https://docs.sympy.org/latest/', None)}
