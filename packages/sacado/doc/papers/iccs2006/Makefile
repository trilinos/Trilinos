paper.dvi: paper.tex paper.bib
	latex paper && bibtex paper && latex paper && latex paper

paper.ps: paper.dvi
	dvips -Ppdf -o paper.ps paper.dvi

paper.pdf: paper.ps
	ps2pdf paper.ps paper.pdf

tex-clean:
	rm -f *.aux *.bbl *.blg *.log *.toc

clean: tex-clean
	rm -f *.ps *.pdf *.dvi
