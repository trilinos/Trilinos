default: esi-spec.ps

esi-spec.ps: esi-spec.dvi
	echo "Creating ps file"
	dvips esi-spec.dvi -o

esi-spec.dvi: *.tex
	echo "Creating dvi file"
	latex esi-spec.tex
	latex esi-spec.tex

clean:
	rm -f *.dvi *.log *.aux *.ps
