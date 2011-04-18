all: KurtisNusbaum.pdf

KurtisNusbaum.pdf: KurtisNusbaum.tex KurtisNusbaum.bib ../OptikaThesis.tex basicTramonto/main.cpp ../OptikaFunctions.tex title.tex ../OptikaIntro.tex
	pdflatex KurtisNusbaum
	bibtex KurtisNusbaum
	pdflatex KurtisNusbaum
	pdflatex KurtisNusbaum

clean:
	rm KurtisNusbaum.pdf KurtisNusbaum.bbl KurtisNusbaum.blg KurtisNusbaum.log KurtisNusbaum.aux *toc
