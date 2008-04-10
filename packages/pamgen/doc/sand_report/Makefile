
# Important note: converting to pdf seems to loose the top space of the
# cover page.  This is due to the tex installation using the European A4 page size by default.
# fix this and the problem will go away.  The -t letter seems to fix this.

REPORT = Pamgen_Report

COMMON_FILES = macros.tex \
               DOEbwlogo.eps snllineblk.eps \
               DOEbwlogo.pdf snllineblk.pdf

LATEX = latex
PDFLATEX = pdflatex
BIBTEX = bibtex
MAKEINDEX = makeindex
DVIPS = dvips -t letter -u pdftex.map
LATEX2HTML = latex2html -split 0 -show_section_numbers -toc_depth 5 -nonavigation -antialias -antialias_text -white

INPUT =     \
pamgen_acknowledgment.tex\
pamgen_introduction.tex\
Pamgen_Report.tex\
pamgen_mesh.tex\
pamgen_rtc.tex\
read_into_memory.tex

SUPPORT_AUX = $(REPORT).aux $(INPUT:.tex=.aux)
SUPPORT_LOG = $(REPORT).log $(INPUT:.tex=.log)

.PHONY: dvi ps pdf all clean realclean cleanall

$(REPORT).dvi: $(REPORT).tex $(INPUT) $(REPORT).bib $(COMMON_FILES)
	$(LATEX) $(REPORT)
	$(LATEX) $(REPORT)
	if test -r $(REPORT).idx ; then \
	  echo "$(MAKEINDEX) $(REPORT)" ; \
	  $(MAKEINDEX) $(REPORT).idx ; \
	fi
	if test -r $(REPORT).bib ; then \
	  echo "$(BIBTEX) $(REPORT)" ; \
	  $(BIBTEX) $(REPORT) ; \
	fi
	$(LATEX) $(REPORT)
	$(LATEX) $(REPORT)

$(REPORT).ps: $(REPORT).dvi
	$(DVIPS) -o $(REPORT).ps $(REPORT)

$(REPORT).pdf: $(REPORT).tex $(INPUT) $(REPORT).bib $(COMMON_FILES)
	$(PDFLATEX) $(REPORT)
	$(PDFLATEX) $(REPORT)
	if test -r $(REPORT).idx ; then \
	  echo "$(MAKEINDEX) $(REPORT).idx" ; \
	  $(MAKEINDEX) $(REPORT).idx ; \
	fi
	if test -r $(REPORT).bib ; then \
	  echo "$(BIBTEX) $(REPORT)" ; \
	  $(BIBTEX) $(REPORT) ; \
	fi
	$(PDFLATEX) $(REPORT)
	$(PDFLATEX) $(REPORT)

html: $(REPORT).tex $(INPUT) $(REPORT).bib $(COMMON_FILES)
	rm -rf $(REPORT)
	$(LATEX) $(REPORT)
	$(LATEX) $(REPORT)
	if test -r $(REPORT).idx ; then \
	  echo "$(MAKEINDEX) $(REPORT).idx" ; \
	  $(MAKEINDEX) $(REPORT).idx ; \
	fi
	if test -r $(REPORT).bib ; then \
	  echo "$(BIBTEX) $(REPORT)" ; \
	  $(BIBTEX) $(REPORT) ; \
	fi
	$(LATEX) $(REPORT)
	$(LATEX) $(REPORT)
	mkdir $(REPORT)
	echo "\batchmode \documentclass{article} \begin{document} \end{document}" > $(REPORT)/images.tex
	$(LATEX2HTML) $(REPORT).tex
	rm -f $(REPORT).aux $(REPORT).toc $(REPORT).idx $(REPORT).ind \
           $(REPORT).ilg $(REPORT).log $(REPORT).out \
           $(REPORT).bbl $(REPORT).blg $(REPORT).lof $(REPORT).lot \
           $(SUPPORT_AUX) $(SUPPORT_LOG)

all: cleanall
	$(MAKE) html
	$(MAKE) $(REPORT).dvi
	$(MAKE) $(REPORT).ps
	$(MAKE) $(REPORT).pdf

dvi: $(REPORT).dvi
	@echo " "
	@echo "Finished"

ps: $(REPORT).ps
	@echo " "
	@echo "Finished"

force:

librarypdf: $(REPORT).ps force
	@echo " "
	ps2pdf $(REPORT).ps
	@echo "Finished"

pdf: $(REPORT).pdf
	@echo " "
	@echo "Finished"

# if there is no bibliography file, this noop target will prevent errors
$(REPORT).bib:
	@echo " " > /dev/null

clean:
	rm -f $(REPORT).aux $(REPORT).toc $(REPORT).idx $(REPORT).ind \
           $(REPORT).ilg $(REPORT).log $(REPORT).out \
           $(REPORT).bbl $(REPORT).blg $(REPORT).lof $(REPORT).lot \
           $(SUPPORT_AUX) $(SUPPORT_LOG)

realclean cleanall: clean
	rm -f $(REPORT).ps $(REPORT).pdf $(REPORT).dvi
	rm -rf $(REPORT)
