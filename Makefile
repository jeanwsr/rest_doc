# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

multilang:
	SPHINX_LANGUAGE=en $(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)"/html/en $(SPHINXOPTS) $(O)
	SPHINX_LANGUAGE=zh_CN $(SPHINXBUILD) -D language=zh_CN -b html "$(SOURCEDIR)_zh" "$(BUILDDIR)"/html/zh_CN $(SPHINXOPTS) $(O)
	cp index.html build/html/index.html
