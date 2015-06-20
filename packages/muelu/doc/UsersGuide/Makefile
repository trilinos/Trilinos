#

#TEX = latex
TEX = pdflatex

BTEX = bibtex

BIBFILE = mueluguide.bib

PAPER = mueluguide

PARTS = \
		SANDdistribution.tex \
		abstract.tex \
		acknowledgment.tex \
		appendix.tex \
		definitions.tex \
		introduction.tex \
		license.tex \
		mueluguide.tex \
		multigrid.tex \
		options.tex \
		options_aggregation.tex \
		options_general.tex \
		options_misc.tex \
		options_multigrid.tex \
		options_rebalancing.tex \
		options_smoothing_and_coarse.tex \
		paramlist.tex \
		paramlist_hidden.tex \
		performance.tex \
		references.tex \
		starting.tex

PICTURES =

all: ${PAPER}.pdf

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
	dvips -Ppdf -o ${PAPER}.ps ${PAPER}

${PAPER}.pdf:    ${PAPER}.dvi

tarfile:
	tar chvf ${PAPER}.tar Makefile ${PARTS} ${PICTURES} ${BIBFILE}
	compress ${PAPER}.tar

clean:
	/bin/rm -f ${PAPER}.pdf ${PAPER}.ps ${PAPER}.dj ${PAPER}.dvi *.aux *.bbl *.blg *.log *.lof *.lot *.toc *.tar *.tar.Z ${PAPER}.out
