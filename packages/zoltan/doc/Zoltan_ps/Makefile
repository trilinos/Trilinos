#
# This Makefile makes PostScript versions of the User's Guide 
# and Developer's Guide from the HTML source.
# 

# Directory for ps files, relative to ug_html and dev_html directories
PSDIR = ../../Zoltan_ps

# Set to your local HTML to PS conversion program
HTML2PS = /Net/proj/zoltan/arch/all/bin/html2ps

# Set the options for your conversion program
HTML2PS_FLAGS = -W b -f $(PSDIR)/Zoltan.html2ps 

all: ug dev

ug:
	(cd ../Zoltan_html/ug_html; \
	$(HTML2PS) $(HTML2PS_FLAGS) -o $(PSDIR)/ug.ps ug.html)

dev:
	(cd ../Zoltan_html/dev_html; \
	$(HTML2PS) $(HTML2PS_FLAGS) -o $(PSDIR)/dev.ps dev.html)
