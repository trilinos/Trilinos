
LATEX = pdflatex
MAKEINDEX = makeindex
MAKEGLOSSARY = makeglossary
BIBTEX = bibtex

DEPS = \
	TempusRequirementsAbstract.tex \
  TempusRequirementsAcknowledgment.tex \
	TempusRequirementsSAND.tex \
	SANDdistribution.tex \
  SANDmath.tex \
  SANDreport.cls \
	

default: help

all: article sand

article: TempusRequirements.pdf

TempusRequirements.pdf: TempusRequirements.tex $(DEPS)
	$(LATEX) TempusRequirements.tex
	$(MAKEINDEX) TempusRequirements
	$(BIBTEX) TempusRequirements
	$(LATEX) TempusRequirements.tex
	$(LATEX) TempusRequirements.tex

sand: TempusRequirementsSAND.pdf

TempusRequirementsSAND.pdf: TempusRequirementsSAND.tex $(DEPS)
	$(LATEX) TempusRequirementsSAND.tex
	$(MAKEINDEX) TempusRequirementsSAND
	$(BIBTEX) TempusRequirementsSAND
	$(LATEX) TempusRequirementsSAND.tex
	$(LATEX) TempusRequirementsSAND.tex

#%.ps: %.dvi
#	dvips -t letter -q $< -o $(<:.dvi=.ps)
#	ps2pdf -dPDFSETTINGS=/prepress $(<:.dvi=.ps)

clean:
	@rm -f \
	*.dvi *.aux *.toc *.idx *.ind *.ilg \
	*.log *.out *.bbl *.blg *.lof *.lot *.tex~ *.log *~

cleanall: clean
	@rm -f TempusRequirements.pdf TempusRequirementsSAND.pdf

help:
	@echo ''
	@echo '  Build TempusRequirements Manuals:'
	@echo ''
	@echo ' Run make with one of the following options:'    
	@echo ''
	@echo ' a) all........build all manuals'
	@echo ' b) clean......remove all built objects'
	@echo ' b) cleanall...clean + remove TempusRequirements.pdf and TempusRequirementsSAND.pdf'
	@echo ' c) article....build manual using article format'
	@echo ' d) sand.......build manual using SAND format'

