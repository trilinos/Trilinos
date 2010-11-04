# ReStructured Text conversion utilities
RST2HTML  = rst2html.py
RST2LATEX = rst2latex.py

HTML_FLAGS = --generator --time --section-numbering --no-xml-declaration
TEX_FLAGS  = --generator --time --section-numbering --use-latex-toc --use-latex-docinfo

# LaTeX
LATEX = pdflatex

# ReStructured Text main source file
SRC_FILE = DevelopersGuide.txt

# Files that $(SRC_FILE) depends on
DEPEND_FILES = AddingModules.txt BuildSystem.txt Documentation.txt \
	       ExceptionHandling.txt Introduction.txt \
	       ReferenceCountedPointers.txt Testing.txt

# Web page to be made
WEB_PAGE = $(SRC_FILE:txt=html)

# PDF file generation
TEX_FILE = $(SRC_FILE:txt=tex)
PDF_FILE = $(TEX_FILE:tex=pdf)

AUX_FILES = $(TEX_FILE:tex=aux) \
	    $(TEX_FILE:tex=out) \
	    $(TEX_FILE:tex=log) \
	    $(TEX_FILE:tex=toc)

# Default target depends on all of the recursive subdirectories
all: $(WEB_PAGE) $(PDF_FILE)

$(WEB_PAGE): $(SRC_FILE) $(DEPEND_FILES)
	$(RST2HTML) $(HTML_FLAGS) $< $@

$(TEX_FILE): $(SRC_FILE) $(DEPEND_FILES)
	$(RST2LATEX) $(TEX_FLAGS) $< $@

$(PDF_FILE): $(TEX_FILE)
	$(LATEX) $<
	$(LATEX) $<
	$(LATEX) $<

clean:
	$(RM) $(TEX_FILE) $(AUX_FILES) *~

clobber: clean
	$(RM) $(WEB_PAGE) $(PDF_FILE)

.PHONY : clean clobber
