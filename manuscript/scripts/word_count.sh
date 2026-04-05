#!/usr/bin/env bash
# Word count for LaTeX manuscript using texcount
set -euo pipefail

if ! command -v texcount &>/dev/null; then
    echo "Error: texcount not found. Install via TeX Live or: brew install texcount"
    exit 1
fi

echo "=== Word Count ==="
texcount -utf8 -inc -total main.tex
echo ""
echo "=== Breakdown by section ==="
texcount -utf8 -inc -sub=section main.tex
