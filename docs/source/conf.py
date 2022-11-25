# Configuration file for the Sphinx documentation builder.

# -- Project information

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

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'nbsphinx',
    'autoapi.extension'
]

autoapi_type = 'python'
autoapi_dirs = [join(SRC_PATH, 'freeflux')]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
