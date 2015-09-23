#

TEX = latex

BTEX = bibtex

BIBFILE = AmesosReferenceGuide.bib

PAPER = AmesosReferenceGuide

PARTS = AmesosReferenceGuide.tex 

PICTURES =

all: ${PAPER}.ps

########################################################################
#
# Dependencies - 3 are required for latex/bibtex
# - a rule must exist for the .dvi, the .bbl and the .aux file
#   in order for a clean make with a new .tex file
#
########################################################################

${PAPER}.dvi: ${PAPER}.tex ${PARTS} ${PAPER}.bbl
	${TEX} ${PAPER}
	${TEX} ${PAPER}

${PAPER}.bbl: ${PAPER}.aux
	${BTEX} ${PAPER}

${PAPER}.aux: ${PAPER}.tex ${PARTS} force
	${TEX} ${PAPER}

force:

${PAPER}.ps:    ${PAPER}.dvi
	dvips -o ${PAPER}.ps ${PAPER}

tarfile:
	tar chvf ${PAPER}.tar Makefile ${PARTS} ${PICTURES} ${BIBFILE} ${SIAM}
	compress ${PAPER}.tar

clean:
	/bin/rm -f *.pdf ${PAPER}.ps ${PAPER}.dj ${PAPER}.dvi *.aux *.bbl *.blg *.log *.lof *.lot *.toc *.tar *.tar.Z


