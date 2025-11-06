# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'REST'
copyright = '2025, REST developers'
author = 'REST developers'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
]

templates_path = ['_templates']
exclude_patterns = []

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']

if html_theme == 'pydata_sphinx_theme':

    html_theme_options = {
        # "logo": {
        #     "text": "REST",
        #     "image_light": "../logo/logo-64x64.png",
        #     "image_dark": "../logo/logo-64x64.png",
        # },
        "icon_links": [
            {
                "name": "Gitee",
                "url": "https://gitee.com/restgroup/rest",
                "icon": "fa-brands fa-git",
                "type": "fontawesome",
            },
        ],
        # "use_edit_page_button": True,
        "show_toc_level": 1,
        "secondary_sidebar_items": {
            "**": ["page-toc"],
        },
        "show_prev_next": False,
        # "navbar_align": "left",
    }
