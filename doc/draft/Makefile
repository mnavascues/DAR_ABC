all: main.pdf diff.pdf

main.pdf: main.tex references.bib
	pdflatex main.tex
	biber main
	pdflatex main.tex
	pdflatex main.tex
	
diff.tex: main.tex
	latexdiff previous_version.tex main.tex > diff.tex

diff.pdf: diff.tex references.bib
	pdflatex diff.tex
	biber diff
	pdflatex diff.tex
	pdflatex diff.tex

clean:
	rm -f *.log *.dvi *.aux
	rm -f *.blg *.bbl *.bcf *.run.xml
	rm -f *.log
	rm -f *.out
	rm -f *.bak
	rm -f diff.tex
