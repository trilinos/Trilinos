all: buildGraphs
	pdflatex matmatreport
	pdflatex matmatreport
	bibtex matmatreport
	pdflatex matmatreport

buildGraphs: justtpetra allmats atrans

justtpetra:
	make -C tpetraDevData all
	
allmats:
	make -C allMatsData all

atrans:
	make -C atransData all


clean:
	rm -f *toc *aux *log matmatreport.pdf *bbl *blg *out *lof *lot
	make -C tpetraDevData clean
	make -C allMatsData clean
	make -C atransData clean
