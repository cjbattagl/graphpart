.DEFAULT_GOAL := all

DRAFT_BUILD_DIR = draft

all: draft-build-dir $(DRAFT_BUILD_DIR)/kdd15_workshop_main.pdf

.PHONY: draft-build-dir
draft-build-dir:
	test -d $(DRAFT_BUILD_DIR)/ || mkdir -p $(DRAFT_BUILD_DIR)

.PHONY: $(DRAFT_BUILD_DIR)/kdd15_workshop_main.pdf
$(DRAFT_BUILD_DIR)/kdd15_workshop_main.pdf:
	latexmk -f -pdf -output-directory=$(DRAFT_BUILD_DIR)/ kdd15_workshop_main.tex

clean:
	rm -rf $(DRAFT_BUILD_DIR)/ *~ core *.blg *.toc *.bbl *.out *.aux *.log

# eof
