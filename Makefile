#
# Makefile for the SAND report Example files
# $Id: Makefile,v 1.11 2008/03/01 16:52:02 rolf Exp $
#
FILES	= OptikaSANDReport

all:	pdf

pdf:	$(addsuffix .pdf, $(FILES))

ps:	$(addsuffix .ps, $(FILES))

# Build one example using pdflatex
OptikaSANDReport.pdf: OptikaSANDReport.tex Optika.bib SANDreport.cls
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<

%.pdf:	%.ps
	ps2pdf13 $< $@

%.ps:	%.dvi
	dvips -Ppdf -o $@ $<

%.dvi:	%.tex Opitka.bib SANDreport.cls
	latex $<
	bibtex $(basename $<)
	latex $<
	latex $<

clean:
	@rm -f $(addsuffix .aux, $(FILES) $(PFILES)) $(addsuffix .bbl, $(FILES) $(PFILES))
	@rm -f $(addsuffix .blg, $(FILES) $(PFILES)) $(addsuffix .lof, $(FILES) $(PFILES))
	@rm -f $(addsuffix .log, $(FILES) $(PFILES)) $(addsuffix .lot, $(FILES) $(PFILES))
	@rm -f $(addsuffix .toc, $(FILES) $(PFILES))
	@rm -f texput.log
	@rm -f Mark*.aux
	@rm -f SANDdistribution.aux

realclean:	clean
	@rm -f $(addsuffix .pdf, $(FILES) $(PFILES)) $(addsuffix .ps, $(FILES) $(PFILES))
	@rm -f $(addsuffix .dvi, $(FILES) $(PFILES)) $(addsuffix .out, $(FILES) $(PFILES))
