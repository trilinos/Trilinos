
LATEX=lualatex
BIBTEX=bibtex
TARGET=ThreadPool

default: thread-pool

thread-pool:
	$(LATEX) $(TARGET)

bib:
	$(BIBTEX) $(TARGET)

.PHONY: clean

clean:
	@rm -f $(TARGET)-blx.bib \
	       $(TARGET).loa \
	       $(TARGET).lof \
	       $(TARGET).lot \
	       $(TARGET).toc \
	       $(TARGET).aux \
	       $(TARGET).idx \
	       $(TARGET).out \
	       $(TARGET).auxlock \
	       $(TARGET).bbl \
	       $(TARGET).blg \
	       $(TARGET).log \
	       $(TARGET).pdf \
	       $(TARGET).xwm \
	       $(TARGET)-figure*.pdf \
	       $(TARGET)-figure*.log \
	       $(TARGET)-figure*.dpth \
	       $(TARGET).snm \
	       $(TARGET).nav
