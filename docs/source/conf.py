# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

### Add pyFuRNAce to the path
import os
import sys

sys.path.insert(0, os.path.abspath("../../"))
###

project = "pyFuRNAce"
copyright = "2025, Luca Monari"
author = "Luca Monari"
release = "0.0.7"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    # "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "myst_parser",
    "nbsphinx",
    "numpydoc",
    "sphinx_copybutton",
    "sphinxcontrib.bibtex",
]

# summary generation
autosummary_generate = True  # Turn on autosummary

# bibtex settings
bibtex_bibfiles = ["refs.bib"]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
    ".ipynb": "jupyter_notebook",
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["static"]
