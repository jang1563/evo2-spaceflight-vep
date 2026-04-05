# bioRxiv Manuscript Harness

[![Build LaTeX](https://github.com/jang1563/bioarxiv-harness/actions/workflows/build.yml/badge.svg)](https://github.com/jang1563/bioarxiv-harness/actions/workflows/build.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A production-ready LaTeX template for bioRxiv preprint submissions. Includes automated builds, pre-submission validation, supplementary materials with cross-references, and a clean single-column layout optimized for readability.

## What Does It Look Like?

<p align="center">
  <img src="img/preview_preprint.png" alt="Preprint layout preview" width="600">
</p>

> Single-column layout with "Preprint --- not peer reviewed" header, line numbering, lead author footer, and ORCID badges. Switch to double-spaced journal submission mode with `make mode-submit`.

## Quick Start

```bash
make              # Build main.pdf + supplement.pdf
make biorxiv      # Merge into single manuscript.pdf for submission
make watch        # Continuous rebuild on file save
make wordcount    # Word count via texcount
make check        # Pre-submission validation (8-point checklist)
make mode-submit  # Switch to journal submission mode (double-spaced)
make clean        # Remove all build artifacts
```

## Requirements

| Tool | Purpose | Install (macOS) |
|------|---------|-----------------|
| **TeX Live 2023+** | pdflatex, latexmk, bibtex | `brew install --cask mactex` |
| **Ghostscript** | Merge PDFs (`make biorxiv`) | `brew install ghostscript` |
| **texcount** | Word counting (`make wordcount`) | Included with TeX Live |
| **git-latexdiff** | Tracked-changes PDF (`make diff`) | `brew install git-latexdiff` (optional) |

## Directory Structure

```
main.tex              Main manuscript
supplement.tex        Supplementary materials (cross-refs to main via xr)
references.bib        Shared bibliography
style/bioarxiv.sty    Style package (header, line numbers, formatting)
sections/             Section files (\input'd by main.tex)
figures/              Figure files (.pdf, .png)
tables/               Table files (optional)
scripts/              Build helpers (word count, submission check)
```

## Usage

1. **Copy** this directory for a new manuscript.
2. **Edit** `main.tex` metadata (title, authors, affiliations, `\leadauthor`).
3. **Fill** section files in `sections/`.
4. **Add** figures to `figures/` and references to `references.bib`.
5. **Build**: `make biorxiv`
6. **Upload** `manuscript.pdf` to [bioRxiv](https://www.biorxiv.org/submit-a-manuscript).

## Features

- Clean single-column layout (11pt, 1-inch margins, 1.5x line spacing)
- **Preprint / journal submission toggle** via `\usepackage[submit]{bioarxiv}`
  - Preprint: 1.5x spacing, no line numbers option
  - Submit: double-spaced with mandatory line numbers
- "Preprint --- not peer reviewed" header with running title
- Lead author footer via `\leadauthor{}`
- Line numbering enabled by default (disable: `\usepackage[nolinenumbers]{bioarxiv}`)
- ORCID support via `\orcidlink{}`
- Biology-specific SI units: `\Molar`, `\rpm`, `\gforce`, `\bp`, `\nt`, `\units`, `\CFU`
- Supplementary materials with S-prefixed counters and cross-references to main
- Pre-submission checklist: TODO markers, undefined refs, file size, line numbers
- CI/CD via GitHub Actions — validates compilation on every push

## Pre-Submission Validation

Run `make check` to validate your manuscript before submission:

```
  [PASS] main.pdf exists (476K)
  [PASS] No TODO/FIXME markers found
  [PASS] No undefined references
  [PASS] No undefined citations
  [PASS] Word count: 4521 (within bioRxiv range)
  [PASS] All referenced figures exist
  [PASS] Line numbers enabled
  [PASS] PDF under 40 MB bioRxiv limit
```

## Overleaf

Upload all files. The Makefile won't work on Overleaf, but latexmk picks up `.latexmkrc` automatically. Set `main.tex` as the main document. Build supplement separately by changing the main document to `supplement.tex`.

## bioRxiv Submission Notes

- bioRxiv accepts **PDF only** for LaTeX manuscripts (no .tex upload)
- Maximum file size: **40 MB**
- No formatting requirements — our style choices are for readability
- Select a CC license during submission (CC BY recommended for preprints)
- Preprints typically appear within 72 hours

## Acknowledgments

Inspired by [kourgeorge/arxiv-style](https://github.com/kourgeorge/arxiv-style) and [quantixed/manuscript-templates](https://github.com/quantixed/manuscript-templates).

## License

[MIT](LICENSE)
