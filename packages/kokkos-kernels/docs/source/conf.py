# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.append(os.path.abspath("./_ext"))
sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------
# The master toctree document.
master_doc = "index"

project = 'Kokkos Kernels'
copyright = 'Copyright Contributors to the Kokkos project'
author = 'lots of people'

# The full version, including alpha/beta/rc tags
# release = '4.5.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["myst_parser",
              "sphinx.ext.autodoc",
              "sphinx.ext.viewcode",
              "sphinx.ext.intersphinx",
              "sphinx_copybutton",
              "sphinx_design",
              "cppkokkos"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "default"
pygments_dark_style = "monokai"

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'furo'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['.']

html_css_files = ['hacks.css']

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

myst_heading_anchors = 4

myst_enable_extensions = [
    "dollarmath",
    "amsmath"
    ]
# need to figure out why this does not work
# rst_prolog = """
# .. include:: special.rst
# """

rst_prolog = """
.. role:: cppkokkos(code)
    :language: cppkokkos

.. role:: cpp(code)
   :language: cpp
"""
