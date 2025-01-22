# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Add the project directory to sys.path
sys.path.insert(0, os.path.abspath('../../'))  # Adjust if your structure differs
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Planner for Transit Observations'
copyright = '2024, Michal Steiner'
author = 'Michal Steiner'
release = '0.1.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',   # Auto-generate API docs
    'sphinx.ext.napoleon',  # Support for Google-style and NumPy-style docstrings
    'sphinx.ext.viewcode',   # Include source code links
    'nbsphinx',
    'sphinx.ext.intersphinx',
]

# Mock specific modules
autodoc_mock_imports = [
    'PTO.simulations', 
    'PTO.transits',
]


# Add intersphinx mapping if needed
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}

# Enable better cross-referencing
autodoc_default_options = {
    'members': True,
    'show-inheritance': True,
    'special-members': '__init__'
}

# Optional: Configure nbsphinx
nbsphinx_execute = 'auto'  # 'always' or 'never'
nbsphinx_allow_errors = True  # Set to False in production

templates_path = ['_templates']
exclude_patterns = []

html_theme_options = {
    # 'navigation_depth': 2,
    # 'collapse_navigation': False  # This ensures the sidebar starts expanded
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']