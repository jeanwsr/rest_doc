# RESTdoc repository

## How to contribute

* Write your awesome documentation in Markdown file under `source/user`, such as `dft.md`, for the usage of certain feature. Add a new file if needed.
* Write corresponding Markdown file in `source/contributor`, if needed, for code structure, brief API doc, etc.
* run `make html` to build the doc. (Need sphinx and dependencies)
* run `python -m http.server --directory build/html` to serve and visit `localhost:8000` in browser to preview the pages.