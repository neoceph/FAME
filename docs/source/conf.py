# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from setuptools import setup, find_packages
from unittest.mock import MagicMock

sys.path.insert(0, os.path.abspath('../../src'))  # Path to your src folder


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FAME'
copyright = '2024, Abdullah Al Amin'
author = 'Abdullah Al Amin'
release = '2024.11'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
autoclass_content = "both"

class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()
    
sys.modules[
    'mpi4py'
    ] = Mock()

on_rtd = os.environ.get('READTHEDOCS') == 'True'

install_requires = [
    'other-dependencies',
]

if not on_rtd:
    install_requires.append('mpi4py')

setup(
    name=project,
    version=release,
    packages=find_packages(),
    install_requires=install_requires,
)