# Configuration file for the Sphinx documentation builder.

import sys
from os.path import dirname, join

SRC_PATH = join(dirname(dirname(dirname(__file__))), 'src')
sys.path.insert(0, SRC_PATH)


# ------------------------ General configuration ------------------------
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.intersphinx',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinx.ext.autosummary',
              'nbsphinx',
              'autoapi.extension']

autoapi_type = 'python'
autoapi_dirs = [join(SRC_PATH, 'freeflux')]
autoapi_add_toctree_entry = False

source_suffix = '.rst'
master_doc = 'index'
# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']


# ------------------------ Project information ------------------------
project = 'freeflux'
copyright = '2022-2024, Chao Wu'
author = 'Chao Wu'
from freeflux import __version__ as version
release = version

pygments_style = 'sphinx'


# ------------------------ Options for HTML output ------------------------
mathjax_path = ('https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML')


# ------------------------ Options for LaTeX output ------------------------
# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, 
#  documentclass [howto, manual, or own class]).
latex_documents = [(master_doc, 
                    project+'.tex', 
                    project+' Documentation',
                    author,
                    'manual')]


# ------------------------ Options for manual page output ------------------------
# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, 
              project, 
              project+' Documentation',
              [author], 
              1)]


# ------------------------ Options for Texinfo output ------------------------
# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [(master_doc,
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
                       'pyomo': ('https://pyomo.readthedocs.io/en/stable/', None),
                       'sympy': ('https://docs.sympy.org/latest/', None)}
