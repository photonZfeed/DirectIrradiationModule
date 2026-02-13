# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Multi-Batch Screening Reactor - Direct Irradiation Module Design'
copyright = '2026, Benedikt Wiedemann, Daniel Kowalczyk et al.'
author = 'Benedikt Wiedemann, Daniel Kowalczyk et al.'

# -- Path setup --------------------------------------------------------------
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    # 'sphinx.ext.autosectionlabel',
    'sphinx_rtd_theme',
]

numfig = True

numfig_format = {
    'figure': 'Figure %s',
    'table': 'Table %s',
    'code-block': 'Code %s',
}

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': False,
}

autodoc_mock_imports = ["numpy", "scipy", "matplotlib", "pandas", "numba", "openpyxl", "tqdm", "ICIW_Plots"]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']