
LATEX = pdflatex
MAKEINDEX = makeindex
MAKEGLOSSARY = makeglossary
BIBTEX = bibtex

DEPS = \
  RythmosConvergenceTestExamples.tex \
  RythmosExplicitStepperFormulation.tex \
  RythmosImplicitStepperFormulation.tex \
  RythmosTheoryBDF.tex \
  RythmosTheoryIRK.tex \
  RythmosTheoryIntroduction.tex \
  RythmosUsersManualParameterList.tex \
  RythmosUsersManualParameterList_Raw.tex \
  rab_commands.tex \
  RythmosTheory.bib \
  ../design_document/RythmosDesignSAND.bib

	

default: help

all: article sand

article: Rythmos.pdf

Rythmos.pdf: Rythmos.tex $(DEPS)
	$(LATEX) Rythmos.tex
	$(MAKEINDEX) Rythmos
	$(BIBTEX) Rythmos
	$(LATEX) Rythmos.tex
	$(LATEX) Rythmos.tex

sand: RythmosSAND.pdf

RythmosSAND.pdf: RythmosSAND.tex $(DEPS)
	$(LATEX) RythmosSAND.tex
	$(MAKEINDEX) RythmosSAND
	$(BIBTEX) RythmosSAND
	$(LATEX) RythmosSAND.tex
	$(LATEX) RythmosSAND.tex

#%.ps: %.dvi
#	dvips -t letter -q $< -o $(<:.dvi=.ps)
#	ps2pdf -dPDFSETTINGS=/prepress $(<:.dvi=.ps)

clean:
	@rm -f RythmosSAND.pdf \
	*.dvi *.aux *.toc *.idx *.ind *.ilg \
	*.log *.out *.bbl *.blg *.lof *.lot *.tex~ *.log *~

help:
	@echo ''
	@echo '  Build Rythmos Manuals:'
	@echo ''
	@echo ' Run make with one of the following options:'    
	@echo ''
	@echo ' a) all........build all manuals'
	@echo ' b) clean......remove all built objects'
	@echo ' c) article....build manual using article format'
	@echo ' d) sand.......build manual using SAND format'

