#
# Makefile for the SMemo_PListDoc memo
#

FILES	= SMemo_PListDoc

all:	pdf ps

pdf:	$(addsuffix .pdf, $(FILES))

ps:	$(addsuffix .ps, $(FILES))

%.pdf:	%.tex smemo.cls
	pdflatex $<
	#bibtex $(basename $<)
	pdflatex $<
	pdflatex $<

%.ps:	%.dvi
	dvips -Ppdf -o $@ $<

%.dvi:	%.tex smemo.cls
	latex $<
	#bibtex $(basename $<)
	latex $<
	latex $<

clean:
	@rm -f $(addsuffix .aux, $(FILES)) $(addsuffix .bbl, $(FILES))
	@rm -f $(addsuffix .blg, $(FILES)) $(addsuffix .lof, $(FILES))
	@rm -f $(addsuffix .log, $(FILES)) $(addsuffix .lot, $(FILES))
	@rm -f $(addsuffix .toc, $(FILES))
	@rm -f texput.log

realclean: clean
	@rm -f $(addsuffix .pdf, $(FILES)) $(addsuffix .ps , $(FILES))
	@rm -f $(addsuffix .dvi, $(FILES)) $(addsuffix .out, $(FILES))
