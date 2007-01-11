# $Id$
#
TEXFILES=apibasics.tex contents.tex introduction.tex examples.tex methods.tex

LATEX=latex
PDFLATEX=pdflatex
LATEX2HTML=latex2html

all:	contents.dvi contents.html contents.pdf

contents.html: $(TEXFILES)
	$(LATEX2HTML) contents.tex

contents.pdf: $(TEXFILES)
	$(PDFLATEX) contents
	$(PDFLATEX) contents

contents.dvi: $(TEXFILES)
	$(LATEX) contents
	$(LATEX) contents

DVICLEAN=contents.toc contents.log contents.aux contents.dvi
PDFCLEAN=contents.pdf
HTMLCLEAN=/bin/rm -rf contents

clean:
	/bin/rm -f $(DVICLEAN) $(PDFCLEAN)
	$(HTMLCLEAN)
