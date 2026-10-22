# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

sys.path.insert(0, os.path.abspath("../.."))
sys.path.insert(0, os.path.abspath("../../ligandqc"))
import importlib.metadata
import pathlib
import tomllib
import ligandqc.const

source_location = pathlib.Path(__file__).parent.parent.parent
if (source_location / "pyproject.toml").exists():
    with open(source_location / "pyproject.toml", "rb") as f:
        __version__ = tomllib.load(f)["project"]["version"]
else:
    __version__ = importlib.metadata.version("package")


project = "ligandqc"
copyright = "2026, Christian Bergeest, Mercedes Alfonso-Prieto, Davide Mandelli, Andrea Rizzi, Giulia Rossetti"
author = "Christian Bergeest, Mercedes Alfonso-Prieto, Davide Mandelli, Andrea Rizzi, Giulia Rossetti"
release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinxcontrib.mermaid",
    "sphinx_click",
    "sphinx.ext.intersphinx",
    # "sphinxcontrib.youtube",
    "sphinxemoji.sphinxemoji",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
    "sphinx_copybutton",
    "sphinx.ext.extlinks",
    "sphinx_design",
    "sphinx_substitution_extensions",
]


templates_path = ["_templates"]
bibtex_bibfiles = ["_static/references.bib"]
bibtex_reference_style = "super"
bibtex_default_style = 'unsrt'
exclude_patterns = []

source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}

napoleon_google_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_static_path = ["_static"]
html_title = "LigandQC Manual"
html_logo = "_static/logo2.png"


html_css_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/fontawesome.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/solid.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/brands.min.css",
]

html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#006f45",
        "color-brand-content": "#006f45",
    },
    "dark_css_variables": {
        "color-brand-primary": "#c57600",
        "color-brand-content": "#c57600",
    },
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/brgstfzj/ligandqc",
            "html": "",
            "class": "fa-brands fa-solid fa-github fa-2x",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/ligandqc/",
            "html": "",
            "class": "fa-solid fa-box fa-2x",
        },
    ],
}

rst_epilog = f"""

.. |THRESH_RSR_GOOD| replace:: {ligandqc.const.THRESH_RSR_GOOD}
.. |THRESH_RSR_BAD| replace:: {ligandqc.const.THRESH_RSR_BAD}
.. |THRESH_RSCC_GOOD| replace:: {ligandqc.const.THRESH_RSCC_GOOD}
.. |THRESH_RSCC_BAD| replace:: {ligandqc.const.THRESH_RSCC_BAD}
.. |THRESH_AVG_OCCUPANCY_GOOD| replace:: {ligandqc.const.THRESH_AVG_OCCUPANCY_GOOD}
.. |THRESH_AVG_OCCUPANCY_BAD| replace:: {ligandqc.const.THRESH_AVG_OCCUPANCY_BAD}
.. |THRESH_ATOMIC_CLASHES_GOOD| replace:: {ligandqc.const.THRESH_ATOMIC_CLASHES_GOOD}
.. |THRESH_ATOMIC_CLASHES_BAD| replace:: {ligandqc.const.THRESH_ATOMIC_CLASHES_BAD}
.. |THRESH_STEREO_CHEMICAL_ERRORS_GOOD| replace:: {ligandqc.const.THRESH_STEREO_CHEMICAL_ERRORS_GOOD}
.. |THRESH_STEREO_CHEMICAL_ERRORS_BAD| replace:: {ligandqc.const.THRESH_STEREO_CHEMICAL_ERRORS_BAD}
.. |THRESH_OUTLIERS_BOND_LENGTH_GOOD| replace:: {ligandqc.const.THRESH_OUTLIERS_BOND_LENGTH_GOOD*100:.0f}%
.. |THRESH_OUTLIERS_BOND_LENGTH_BAD| replace:: {ligandqc.const.THRESH_OUTLIERS_BOND_LENGTH_BAD*100:.0f}%
.. |THRESH_OUTLIERS_BOND_ANGLE_GOOD| replace:: {ligandqc.const.THRESH_OUTLIERS_BOND_ANGLE_GOOD*100:.0f}%
.. |THRESH_OUTLIERS_BOND_ANGLE_BAD| replace:: {ligandqc.const.THRESH_OUTLIERS_BOND_ANGLE_BAD*100:.0f}%
.. |THRESH_RMSZ_BOND_LENGTH_BAD| replace:: {ligandqc.const.THRESH_RMSZ_BAD_UPPER}
.. |THRESH_RMSZ_BOND_ANGLE_BAD| replace:: {ligandqc.const.THRESH_RMSZ_BAD_UPPER}
.. |CATEGORY_GOOD| replace:: {ligandqc.const.CATEGORY_GOOD}
.. |CATEGORY_OK| replace:: {ligandqc.const.CATEGORY_OK}
.. |CATEGORY_BAD| replace:: {ligandqc.const.CATEGORY_BAD}
.. |CATEGORY_UNKNOWN| replace:: {ligandqc.const.CATEGORY_UNKNOWN}

"""