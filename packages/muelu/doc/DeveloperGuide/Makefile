MAIN = devGuide

DEPS = $(MAIN) abstract introduction definitions namingConventions codingConventions testing hints commonErrors development gettingStarted

LATEX  = pdflatex
BIBTEX = bibtex

all:	$(MAIN).pdf

$(MAIN).pdf: $(addsuffix .tex, $(DEPS) ) $(addsuffix .bib, $(MAIN) )
	$(LATEX) $<
	$(BIBTEX) $(basename $<)
	$(LATEX) $<
	$(LATEX) $<

clean:
	@rm -f $(addsuffix .aux, $(MAIN) ) $(addsuffix .bbl, $(MAIN) )
	@rm -f $(addsuffix .blg, $(MAIN) ) $(addsuffix .lof, $(MAIN) )
	@rm -f $(addsuffix .log, $(MAIN) ) $(addsuffix .lot, $(MAIN) )
	@rm -f $(addsuffix .toc, $(MAIN) ) $(addsuffix .out, $(MAIN) )
	@rm -f $(addsuffix .pdf, $(MAIN) )
