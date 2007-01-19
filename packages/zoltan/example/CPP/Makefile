# $Id$
#

EXAMPLES = zCPPExample1 zCPPExample2 simpleRCB 

ifndef ZOLTAN_LIB_DIR
all:
	@cd .. && make CPP_Examples
else
all:    $(EXAMPLES)
endif


% : %.cpp ../lib/libexzoltan.a $(ZOLTAN_LIB_DIR)/libzoltan.a
	@echo "Building $@"
	@$(CC) $(CFLAGS) -o $@ $(INCLUDE) $(LIBDIR) $< $(LIBS)

simpleRCB: rectangularMesh.h simpleRCB.cpp
	@echo "Building simpleRCB"
	@$(CC) $(CFLAGS) -o simpleRCB $(INCLUDE) $(LIBDIR) simpleRCB.cpp $(LIBS)

../lib/libexzoltan.a:
	@cd ..; $(MAKE) exampleLibrary

clean:
	rm -f *.o $(EXAMPLES)



