# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os

# Detect which language we're building
current_language = os.environ.get('SPHINX_LANGUAGE', 'en')

if current_language == 'zh_CN':
    source_dir = 'source_zh'
    language = 'zh_CN'
    # html_title = "中文文档"
else:
    source_dir = 'source_en' 
    language = 'en'
    # html_title = "My Documentation"

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

# language = 'en'
locale_dirs = ['locale/']
html_context = {
    'languages': [
        ('en', 'English'),
        ('zh_CN', '中文'),
    ]
}
# html_context = {
#   'current_version' : "1.0",
#   'versions' : [["1.0", "link to 1.0"], ["2.0", "link to 2.0"]],
#   'current_language': 'en',
#   'languages': [["en", "link to en"], ["de", "link to de"]]
# }


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']
html_css_files = [
    'language-switcher.css',
]
html_js_files = [
    'language-switcher.js',
]

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
