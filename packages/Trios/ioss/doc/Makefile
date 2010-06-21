.SUFFIXES: .tex .pdf

all: IOSystem.pdf cth_presto_one_way_coupling.pdf

clean:
	rm -f *.aux *.log *.toc *.aux *.dvi *.out *.pdf *~

IOSystem.pdf: Block_Diagram_Model.png
IOSystem.pdf: DatabaseIO__inherit__graph.png
IOSystem.pdf: GroupingEntity__coll__graph.png
IOSystem.pdf: GroupingEntity__inherit__graph.png
IOSystem.pdf: Ioss_Region.png
IOSystem.pdf: VariableType__inherit__graph.png

.tex.pdf:
	pdflatex $*.tex
	pdflatex $*.tex
