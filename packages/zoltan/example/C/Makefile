# $Id$
#

CEXAMPLES = zoltanSimple zoltanExample1 

all:    $(CEXAMPLES)

% : %.c ../lib/libexzoltan.a 
	$(CC) -o $@ $(INCLUDE) $(LIBDIR) $< $(LIBS)

../lib/libexzoltan.a:
	cd ..; $(MAKE) exampleLibrary

clean:
	rm -f *.o $(CEXAMPLES)



