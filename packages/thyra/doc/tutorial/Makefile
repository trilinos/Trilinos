
.SECONDARY. : UsersGuide.dvi

UsersGuide.pdf : UsersGuide.ps
	ps2pdf  $<

%.ps : %.dvi
	dvips -o $@ $<

%.dvi : %.tex
	latex $<

clean:
	/bin/rm *.aux
	/bin/rm *.log
	/bin/rm *.dvi
	/bin/rm *.ps
	/bin/rm *.pdf

dvi:
	latex UsersGuide
	
	
