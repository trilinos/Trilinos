SHELL = /bin/sh

setup:
	@(if [ ! -f ml_defs.h ]; \
        then \
	   cp ml_defs.H ml_defs.h; \
	   cp ml_include.H ml_include.h; \
        fi;) 

clean:
	@echo "junk" > dummy.h
	@rm -f *.h 

