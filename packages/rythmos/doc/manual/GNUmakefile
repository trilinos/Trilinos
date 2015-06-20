
LATEX=lualatex
BIBTEX=bibtex
TARGET0=Rythmos
TARGET1=RythmosSAND

default: rythmos rythmos-sand

rythmos:
	$(LATEX) $(TARGET0)

rythmos-sand:
	$(LATEX) $(TARGET1)

bib:
	$(BIBTEX) $(TARGET0)

.PHONY: clean

clean:
	@rm -f $(TARGET0)-blx.bib \
	       $(TARGET0).lof \
	       $(TARGET0).lot \
	       $(TARGET0).toc \
	       $(TARGET0).aux \
	       $(TARGET0).idx \
	       $(TARGET0).out \
	       $(TARGET0).auxlock \
	       $(TARGET0).bbl \
	       $(TARGET0).blg \
	       $(TARGET0).log \
	       $(TARGET0).pdf \
	       $(TARGET0).xwm \
	       $(TARGET0)-figure*.pdf \
	       $(TARGET0)-figure*.log \
	       $(TARGET0)-figure*.dpth \
	       $(TARGET0).snm \
	       $(TARGET0).nav
	@rm -f $(TARGET1)-blx.bib \
	       $(TARGET1).lof \
	       $(TARGET1).lot \
	       $(TARGET1).toc \
	       $(TARGET1).aux \
	       $(TARGET1).idx \
	       $(TARGET1).out \
	       $(TARGET1).auxlock \
	       $(TARGET1).bbl \
	       $(TARGET1).blg \
	       $(TARGET1).log \
	       $(TARGET1).pdf \
	       $(TARGET1).xwm \
	       $(TARGET1)-figure*.pdf \
	       $(TARGET1)-figure*.log \
	       $(TARGET1)-figure*.dpth \
	       $(TARGET1).snm \
	       $(TARGET1).nav
	@rm -f x.log
