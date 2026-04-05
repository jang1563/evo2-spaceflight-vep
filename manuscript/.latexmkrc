$pdflatex = 'pdflatex -interaction=nonstopmode -synctex=1 %O %S';
$pdf_mode = 1;
$bibtex_use = 2;
@default_files = ('main.tex');

# Allow LaTeX to find style/bioarxiv.sty
ensure_path('TEXINPUTS', './style//');
