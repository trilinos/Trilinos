# $Id$
#

EXAMPLES = zCPPExample1 zCPPExample2

all:    $(EXAMPLES)

% : %.cpp ../lib/libexzoltan.a 
	$(CC) -o $@ $(INCLUDE) $(LIBDIR) $< $(LIBS)

../lib/libexzoltan.a:
	cd ..; $(MAKE) exampleLibrary

clean:
	rm -f *.o $(EXAMPLES)



