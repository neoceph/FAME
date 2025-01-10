# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath('../../fame'))  # Path to src folder
sys.path.insert(0, os.path.abspath('../..'))  # Path to root

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FAME'
copyright = '2024~2025, Abdullah Al Amin'
author = 'Abdullah Al Amin'
release = '2025.01'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
    'sphinx.ext.mathjax',
    'sphinxcontrib.cairosvgconverter',
]

# Image conversion
# Specify the converter for handling .svg files
image_converter = 'cairosvg'

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = [
    'style.css',  # Specify your custom CSS file
]
autoclass_content = "both"


# conf.py

# Mock heavy or unnecessary imports
autodoc_mock_imports = ["petsc4py", "petsc4py", "jax"]
