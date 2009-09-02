#

TEX = pdflatex

BTEX = bibtex

BIBFILE = 

PAPER = inout

PARTS = inout.tex 

PICTURES =

all: ${PAPER}.pdf

${PAPER}.pdf: ${PAPER}.tex ${PARTS} ${PAPER}.bbl
	${TEX} ${PAPER}

${PAPER}.bbl: ${PAPER}.aux
	${BTEX} ${PAPER}

${PAPER}.aux: ${PAPER}.tex ${PARTS} 
	${TEX} ${PAPER}

clean:
	/bin/rm -f ${PAPER}.ps ${PAPER}.dj ${PAPER}.dvi *.aux *.bbl *.blg *.log *.lof *.lot *.toc *.tar *.tar.Z


