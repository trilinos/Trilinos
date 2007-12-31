# ReStructured Text
RST2HTML = rst2html.py
RFLAGS = --generator --time --no-xml-declaration

# Web pages that need to be made
WEB_PAGES = AddingModules.html BuildSystem.html Documentation.html \
	    ExceptionHandling.html Testing.html

# List all of the subdirectories here for recursive make
SUBDIRS = 

# Default target depends on all of the recursive subdirectories
all: $(WEB_PAGES) $(SUBDIRS)

%.html: %.txt
	$(RST2HTML) $(RFLAGS) $< $@

clean:
	@for dir in $(SUBDIRS); do  \
	    echo Entering $$dir;    \
	    cd $$dir && make clean; \
	    echo Leaving $$dir;     \
	    cd ..;                  \
	 done

# Each recursive make directory needs its own target
#development:
#	cd $@ && $(MAKE)

# The subdirectories are not actually "made", but they do exist, so
# they need to be declared phony targets
.PHONY : $(SUBDIRS) clean
