# RESTdoc repository

## How to contribute

* Write your awesome documentation in Markdown file under `source/user`, such as `dft.md`, for the usage of certain feature. Add a new file if needed.
* Write chinese version in `source_zh`. Write in only one language is ok, just left another language version as "to be added".
* Write corresponding Markdown file in `source/contributor`, if needed, for code structure, brief API doc, etc.
* Keep the Chinese/English wording consistent with the terminology table at `references/terminology.json` (format `"中文": "English"`); add new terms there.
* run `make multilang` to build the doc. (Need sphinx and dependencies, run `pip install -r requirements.txt` if not installed)
* run `python -m http.server --directory build/html` to serve and visit `localhost:8000` in browser to preview the pages.
