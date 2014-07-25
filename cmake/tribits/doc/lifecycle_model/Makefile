
BASE_NAME = TribitsLifecycleModel
LATEX_COMMAND = pdflatex
MAKEINDEX_COMMAND = makeindex
MAKEGLOSSARY_COMMAND = makeglossary
BIBTEX_COMMAND = bibtex

SOURCEFILES= \
  $(BASE_NAME).tex \
  references.bib

#
# Shortcut targets
#

default: $(BASE_NAME).pdf

clean:
	rm -f *.ps $(BASE_NAME).pdf *.dvi *.aux *.toc *.idx *.ind *.ilg \
	*.log *.out *.bbl *.blg *.lof *.lot *.tex~ .log *~

#
# File targets
#

$(BASE_NAME).pdf: $(SOURCEFILES)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(BIBTEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)
