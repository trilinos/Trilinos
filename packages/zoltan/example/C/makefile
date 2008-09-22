# $Id$
#
# Examples, written in C, of parallel applications that use the Zoltan library
# to partition a collection of objects. 
#
# simpleBLOCK.c - The most basic example, with 25 objects, each of which
#   has a global ID and a weight.  The partitioning method is BLOCK, a
#   method which exists as a guide to Zoltan developers on how to write a
#   Zoltan partitioning method.  BLOCK is not expected to be useful in
#   most real world problems.
#
# simpleRCB.c - Recursive coordinate bisection, a geometric method, with
#   25 objects, each of which has a global ID, a weight, and geometric
#   coordinates.
#
# simpleGRAPH.c - ParMETIS, graph partitioning, a graph of 25 vertices
#   with global IDs, connectivity, and vertex weights.
#
# simpleHG.c - Hypergraph partitioning, a hypergraph with 25 hyperedges,
#   each of which has a globalID, an edge weight, and a list of vertices
#   in the hyperedge.  The vertices are the objects to be partitioned,
#   with the goal of achieving a load balance with few hyperedge cuts.
#   Each vertex has a global ID and connectivity.  Vertex weights are
#   supported, but we don't use them in this example.
#
# zoltanExample1 - Partition a simple mesh using three different
#   geometric methods (RCB, recursive inertial bisection, and
#   Hilbert space filling curve) and use Zoltan_Timer routines to
#   time the computation.
#
# zoltanSimple - The first example we wrote, RCB partitioning of
#   a rectangular mesh.
#
# phgExample - A hypergraph example which is really more of a test of
#   hypgergraph parameters and query functions than it is an example.  
#   The hypergraph is less trivial than the 25 edge graph, it has 
#   isolated vertices (vertices not in any hyperedge), some dense edges, 
#   and it's representation forms a picture which makes it easier to 
#   locate errors when debugging new Zoltan hypergraph features.
#

#
# Run these however you run MPI applications:  mpirun -np 4 simpleRCB
#

ifdef PARMETIS_LIBPATH
GRAPH_EXAMPLE=simpleGRAPH
endif

CEXAMPLES = simpleBLOCK simpleRCB $(GRAPH_EXAMPLE) simpleHG zoltanSimple zoltanExample1 subDirs phgExample


LIBZOLTAN = ../lib/libexzoltan.a $(ZOLTAN_LIB_DIR)/libzoltan.a

DOTDOTINCLUDE= $(subst -I../, -I../../, $(INCLUDE))
DOTDOTLIBDIR= $(subst -L../, -L../../, $(LIBDIR))

ifndef ZOLTAN_LIB_DIR
all:
	@cd .. && make C_Examples
else
all:    $(CEXAMPLES)
endif

simpleBLOCK : simpleBLOCK.c simpleGraph.h simpleQueries.h $(LIBZOLTAN)
	@echo "Building $@"
	@$(CC) $(CFLAGS) -o $@ $(INCLUDE) $(LIBDIR) $< $(LIBS)

simpleHG : simpleHG.c simpleGraph.h simpleQueries.h $(LIBZOLTAN)
	@echo "Building $@"
	@$(CC) $(CFLAGS) -o $@ $(INCLUDE) $(LIBDIR) $< $(LIBS)

simpleRCB : simpleRCB.c simpleGraph.h simpleQueries.h $(LIBZOLTAN)
	@echo "Building $@"
	@$(CC) $(CFLAGS) -o $@ $(INCLUDE) $(LIBDIR) $< $(LIBS)

simpleGRAPH : simpleGRAPH.c simpleGraph.h simpleQueries.h $(LIBZOLTAN)
	@echo "Building $@"
	@$(CC) $(CFLAGS) -o $@ $(INCLUDE) $(LIBDIR) $< $(LIBS)

% : %.c $(LIBZOLTAN)
	@echo "Building $@"
	@$(CC) $(CFLAGS) -o $@ $(INCLUDE) $(LIBDIR) $< $(LIBS)

../lib/libexzoltan.a:
	@cd ..; $(MAKE) exampleLibrary

subDirs:
	@cd sparse_matrix; \
        $(MAKE) all "CC=$(CC)" "CFLAGS=$(CFLAGS)" \
                    "INCLUDE=$(DOTDOTINCLUDE)" \
                    "LIBDIR=$(DOTDOTLIBDIR)" "LIBS=$(LIBS)"
clean:
	rm -f *.o  $(CEXAMPLES)
	@cd sparse_matrix; $(MAKE) clean



