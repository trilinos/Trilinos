TEXFILES = build-precond.tex

all: build-precond.pdf

build-precond.pdf: $(TEXFILES) build-precond.bbl
	pdflatex build-precond.tex

build-precond.bbl: build-precond.bib
	make clean
	pdflatex build-precond.tex
	bibtex build-precond
	pdflatex build-precond.tex
	
clean:
	-rm *.blg *.aux *.pdf *.bbl *.log *.out
