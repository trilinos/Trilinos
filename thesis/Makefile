files = KurtisNusbaum.tex KurtisNusbaum.bib ../OptikaThesis.tex basicTramonto/main.cpp ../OptikaFunctions.tex title.tex ../OptikaIntro.tex \
../OptikaAbstract.tex basicXMLTramonto/main.cpp basicXMLTramonto/inputs.xml \
basicTramonto/main.cpp


all: KurtisNusbaum.pdf

KurtisNusbaum.pdf: $(files)
	pdflatex KurtisNusbaum
	bibtex KurtisNusbaum
	pdflatex KurtisNusbaum
	pdflatex KurtisNusbaum

clean:
	rm KurtisNusbaum.pdf KurtisNusbaum.bbl KurtisNusbaum.blg KurtisNusbaum.log KurtisNusbaum.aux *toc
