MAIN   = main
SUPP   = supplement
SHELL  = /bin/bash

SECTIONS = $(wildcard sections/*.tex)
FIGURES  = $(wildcard figures/*.pdf) $(wildcard figures/*.png)
DEPS     = $(MAIN).tex references.bib $(SECTIONS) $(FIGURES)

# ── Build ────────────────────────────────────────────────────────────────────

.PHONY: all main supplement biorxiv watch wordcount check diff clean distclean
.PHONY: mode-preprint mode-submit

all: main supplement

main: $(DEPS)
	latexmk $(MAIN).tex

# supplement depends on main.aux for xr cross-references
supplement: main
	latexmk $(SUPP).tex

# ── bioRxiv Submission ───────────────────────────────────────────────────────

biorxiv: all
	gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pdfwrite \
		-sOutputFile=manuscript.pdf $(MAIN).pdf $(SUPP).pdf
	@echo "Created manuscript.pdf (main + supplement merged)"

# ── Live Preview ─────────────────────────────────────────────────────────────

watch:
	latexmk -pvc $(MAIN).tex

# ── Word Count ───────────────────────────────────────────────────────────────

wordcount:
	@bash scripts/word_count.sh

# ── Pre-Submission Check ─────────────────────────────────────────────────────

check: all
	@bash scripts/check_submission.sh

# ── Revision Diff ────────────────────────────────────────────────────────────

diff:
	@command -v git-latexdiff >/dev/null 2>&1 || \
		{ echo "git-latexdiff not installed. Install: brew install git-latexdiff"; exit 1; }
	git latexdiff HEAD~1

# ── Switch Layout Mode ───────────────────────────────────────────────────────

mode-preprint:
	@sed -i '' 's/\\usepackage\[submit\]{bioarxiv}/\\usepackage{bioarxiv}/' $(MAIN).tex
	@echo "Switched to preprint mode (single-column, 1.5x spacing)"

mode-submit:
	@sed -i '' 's/\\usepackage{bioarxiv}/\\usepackage[submit]{bioarxiv}/' $(MAIN).tex
	@echo "Switched to journal submission mode (double-spaced, line numbers)"

# ── Clean ────────────────────────────────────────────────────────────────────

clean:
	latexmk -C
	latexmk -C $(SUPP).tex
	rm -f manuscript.pdf
	@echo "Cleaned all build artifacts"

distclean: clean
	rm -f $(MAIN).pdf $(SUPP).pdf
