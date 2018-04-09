
PDFENGINE=pdflatex
BIBENGINE=bibtex

TEXFILES=\
	ShyLU_Intro_Tacho.tex \
	ShyLU_Intro_Basker.tex \
	ShyLU_Intro_FastILU.tex \
	ShyLU_Intro_HTS.tex \
	ShyLU_Intro_FROSch.tex \
	ShyLU_Intro_BDDC.tex \
	ShyLU_Intro_Introduction.tex \
	ShyLU_Intro_Ack.tex \
	ShyLU_Intro.tex \
	ShyLU_Intro.bib

ShyLU_Intro.pdf: $(TEXFILES)
	$(PDFENGINE) ShyLU_Intro.tex
	$(BIBENGINE) ShyLU_Intro	
	$(PDFENGINE) ShyLU_Intro.tex
	$(BIBENGINE) ShyLU_Intro	
	$(PDFENGINE) ShyLU_Intro.tex

clean:
	rm *.pdf *.aux *.log *.toc *.out
