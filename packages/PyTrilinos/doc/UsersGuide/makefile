
LATEX=lualatex
BIBTEX=bibtex
TARGET=guide

default: guide

guide:
	$(LATEX) $(TARGET)

bib:
	$(BIBTEX) $(TARGET)

clean:
	@rm -f $(TARGET)-blx.bib \
	       $(TARGET).aux \
	       $(TARGET).toc \
	       $(TARGET).out \
	       $(TARGET).auxlock \
	       $(TARGET).bbl \
	       $(TARGET).blg \
	       $(TARGET).log \
	       $(TARGET).pdf \
	       $(TARGET).xwm \
	       $(TARGET).snm \
	       $(TARGET).nav
	@rm -f figures/*-converted-to.pdf
	@rm -f missfont.log \
	       x.log
