
BASE_NAME = TrilinosSoftwareEngineeringImprovements
LATEX_COMMAND = latex
BIBTEX_COMMAND = bibtex

SOURCEFILES= \
            $(BASE_NAME).tex \
            references.bib

#
# Shortcut targets
#

default: pdf

dvi: $(BASE_NAME).dvi

ps: $(BASE_NAME).ps

pdf: $(BASE_NAME).pdf

clean:
	rm -f *.ps $(BASE_NAME).pdf *.dvi *.aux *.toc *.idx *.ind *.ilg \
	*.log *.out *.bbl *.blg *.lof *.lot *.tex~ .log *~

#
# File targets
#

$(BASE_NAME).dvi: $(SOURCEFILES)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(BIBTEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)
	$(LATEX_COMMAND) $(BASE_NAME)

$(BASE_NAME).ps: $(BASE_NAME).dvi
	dvips -t letter -o $(BASE_NAME).ps $(BASE_NAME)

$(BASE_NAME).pdf: $(BASE_NAME).ps
	ps2pdf $(BASE_NAME).ps $(BASE_NAME).pdf

