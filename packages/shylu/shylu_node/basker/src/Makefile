#!
# Makefile for Baske
# Written by Joshua Dennis Booth
# Sandia National Labs
# Spring 2015
#----------------------------------------------------

include ../make.inc

OBJ = basker.o

all:  basker_lib

basker_lib: $(OBJ)
	$(ARCH) $(ARCHFLAGS) ../lib/libbasker.a $(OBJ)
	$(RANLIB) ../lib/libbasker.a

basker.o: basker.cpp
	$(CC) -g -c $(DFLAGS) $(OPT) $(OMPINC) $(CCINC) $(CCLD) $<

clean:
	rm -f *.o ../lib/libbasker.a


