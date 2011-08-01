REPORT_NAME= matmatreport

all: buildGraphs
	pdflatex $(REPORT_NAME)
	bibtex $(REPORT_NAME)
	pdflatex $(REPORT_NAME)
	pdflatex $(REPORT_NAME)

buildGraphs: justtpetra allmats atrans

justtpetra:
	make -C tpetraDevData all
	
allmats:
	make -C allMatsData all

atrans:
	make -C atransData all


clean:
	rm -f $(REPORT_NAME).toc \
	$(REPORT_NAME).aux \
	$(REPORT_NAME).log \
	$(REPORT_NAME).pdf \
	$(REPORT_NAME).bbl \
	$(REPORT_NAME).blg \
	$(REPORT_NAME).out \
	$(REPORT_NAME).lof \
	$(REPORT_NAME).lot
	make -C tpetraDevData clean
	make -C allMatsData clean
	make -C atransData clean
