import os
import sys
sys.path.insert(0, os.path.abspath('../../drugintfinder'))

# -- General configuration ---------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx_click.ext',
]

templates_path = ['_templates']

source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'Druggable Interactor Finder'
copyright = "2021, Bruce Schultz"
author = "Bruce Schultz"

# The version info for the project you're documenting, acts as replacement
# for |version| and |release|, also used in various other places throughout
# the built documents.
#
# The short X.Y version.
version = '0.1'

# The full version, including alpha/beta/rc tags.
release = '0.2.3'

language = 'en'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'sqlalchemy': ('https://docs.sqlalchemy.org/en/14/', None),
}

# -- Options for HTML output -------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_static_path = []
htmlhelp_basename = 'drugintfinderdoc'


# -- Options for LaTeX output ------------------------------------------

latex_elements = {}
latex_documents = [
    (master_doc, 'drugintfinder.tex',
     'Druggable Interactor Finder Documentation',
     'Bruce Schultz', 'manual'),
]


# -- Options for manual page output ------------------------------------
man_pages = [
    (master_doc, 'drugintfinder',
     'Druggable Interactor Finder Documentation',
     [author], 1)
]

texinfo_documents = [
    (master_doc, 'drugintfinder',
     'Druggable Interactor Finder Documentation',
     author,
     'drugintfinder',
     'One line description of project.',
     'Miscellaneous'),
]
