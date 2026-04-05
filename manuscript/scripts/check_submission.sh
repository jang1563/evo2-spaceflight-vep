#!/usr/bin/env bash
# Pre-submission validation checklist for bioRxiv
set -euo pipefail

PASS=0
FAIL=0
WARN=0

pass() { echo "  [PASS] $1"; ((PASS++)); }
fail() { echo "  [FAIL] $1"; ((FAIL++)); }
warn() { echo "  [WARN] $1"; ((WARN++)); }

echo "=== bioRxiv Submission Checklist ==="
echo ""

# 1. PDF exists and is non-empty
if [ -s main.pdf ]; then
    pass "main.pdf exists ($(du -h main.pdf | cut -f1))"
else
    fail "main.pdf missing or empty — run 'make' first"
fi

# 2. No TODO markers in .tex files
TODO_COUNT=$(grep -r -c 'TODO' sections/ main.tex supplement.tex 2>/dev/null | awk -F: '{s+=$2} END {print s+0}')
if [ "$TODO_COUNT" -eq 0 ]; then
    pass "No TODO markers found"
else
    warn "$TODO_COUNT TODO markers remain in .tex files"
fi

# 3. No undefined references
if [ -f main.log ]; then
    UNDEF_REFS=$(grep -c "LaTeX Warning.*Reference.*undefined" main.log 2>/dev/null || true)
    if [ "$UNDEF_REFS" -eq 0 ]; then
        pass "No undefined references"
    else
        fail "$UNDEF_REFS undefined reference(s) — check main.log"
    fi
else
    warn "main.log not found — build first"
fi

# 4. No missing citations
if [ -f main.log ]; then
    UNDEF_CITES=$(grep -c "Citation.*undefined" main.log 2>/dev/null || true)
    if [ "$UNDEF_CITES" -eq 0 ]; then
        pass "No undefined citations"
    else
        fail "$UNDEF_CITES undefined citation(s) — check main.log"
    fi
else
    warn "main.log not found — build first"
fi

# 5. Word count (informational — bioRxiv has no word limit)
if command -v texcount &>/dev/null; then
    WC=$(texcount -utf8 -inc -total -brief main.tex 2>/dev/null | grep -oE '[0-9]+' | head -1)
    pass "Word count: ~${WC:-unknown} words (bioRxiv has no limit)"
else
    warn "texcount not installed — skipping word count"
fi

# 6. Check that included figures exist
MISSING_FIGS=0
for fig in $(grep -oP '\\includegraphics[^{]*\{[^}]*\}' main.tex sections/*.tex 2>/dev/null | grep -oP '\{[^}]*\}' | tr -d '{}'); do
    # Check with and without figures/ prefix
    if [ ! -f "$fig" ] && [ ! -f "figures/$fig" ]; then
        fail "Missing figure: $fig"
        ((MISSING_FIGS++))
    fi
done
if [ "$MISSING_FIGS" -eq 0 ]; then
    pass "All referenced figures found (or none referenced yet)"
fi

# 7. Line numbers enabled
if grep -q '\\linenumbers' style/bioarxiv.sty 2>/dev/null; then
    if grep -q '\\nolinenumbers' main.tex 2>/dev/null; then
        warn "Line numbers disabled in main.tex (\\nolinenumbers found)"
    else
        pass "Line numbers enabled"
    fi
else
    fail "\\linenumbers not found in style/bioarxiv.sty"
fi

# 8. PDF file size < 40 MB (bioRxiv limit)
if [ -f main.pdf ]; then
    SIZE_BYTES=$(wc -c < main.pdf | tr -d ' ')
    SIZE_MB=$((SIZE_BYTES / 1048576))
    if [ "$SIZE_MB" -lt 40 ]; then
        pass "PDF size: ${SIZE_MB} MB (limit: 40 MB)"
    else
        fail "PDF too large: ${SIZE_MB} MB (bioRxiv limit: 40 MB)"
    fi
fi

# Summary
echo ""
echo "=== Summary: $PASS passed, $FAIL failed, $WARN warnings ==="
if [ "$FAIL" -gt 0 ]; then
    echo "Fix failures before submitting."
    exit 1
fi
