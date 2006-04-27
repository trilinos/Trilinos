#

TEX = pdflatex

BTEX = bibtex

BIBFILE = 

PAPER = MatrixPortal

PARTS = MatrixPortal.tex

PICTURES =

all: MatrixPortal.pdf

MatrixPortal.bbl: MatrixPortal.aux
	${BTEX} MatrixPortal.tex

MatrixPortal.aux: MatrixPortal.tex
	${TEX} MatrixPortal.tex

MatrixPortal.pdf: MatrixPortal.aux
	${TEX} MatrixPortal.tex

clean:
	/bin/rm -f MatrixPortal.pdf *.aux *.bbl *.blg *.log *.lof *.lot *.toc *.tar *.tar.Z


