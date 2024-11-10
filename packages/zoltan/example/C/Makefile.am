
# The following line helps the test harness recover from build errors.

all-local:

include $(top_builddir)/Makefile.export.zoltan

EXEEXT = .exe

noinst_PROGRAMS = simpleBLOCK simpleRCB simpleGRAPH simplePHG problemGRAPH

simpleBLOCK_SOURCES      = simpleBLOCK.c
simpleBLOCK_DEPENDENCIES = $(ZOLTAN_DEPS)
simpleBLOCK_CFLAGS     = $(ZOLTAN_INCLUDES)
simpleBLOCK_LDADD        = $(ZOLTAN_LIBS)

simpleRCB_SOURCES      = simpleRCB.c
simpleRCB_DEPENDENCIES = $(ZOLTAN_DEPS)
simpleRCB_CFLAGS     = $(ZOLTAN_INCLUDES)
simpleRCB_LDADD        = $(ZOLTAN_LIBS)

simpleGRAPH_SOURCES      = simpleGRAPH.c
simpleGRAPH_DEPENDENCIES = $(ZOLTAN_DEPS)
simpleGRAPH_CFLAGS     = $(ZOLTAN_LIBRARY_INCLUDES)
simpleGRAPH_LDADD        = $(ZOLTAN_LIBS)

simplePHG_SOURCES      = simplePHG.c
simplePHG_DEPENDENCIES = $(ZOLTAN_DEPS)
simplePHG_CFLAGS     = $(ZOLTAN_INCLUDES)
simplePHG_LDADD        = $(ZOLTAN_LIBS)

problemGRAPH_SOURCES      = problemGRAPH.c
problemGRAPH_DEPENDENCIES = $(ZOLTAN_DEPS)
problemGRAPH_CFLAGS     = $(ZOLTAN_LIBRARY_INCLUDES)
problemGRAPH_LDADD        = $(ZOLTAN_LIBS)

EXTRA_DIST = classicMakefile objects.txt hypergraph.txt graph.txt mesh.txt

FNAMES=$(EXTRA_DIST:%=$(srcdir)/%)

all-local:
	if test -e Makefile.am ; then : ; else cp $(FNAMES) $(builddir); fi
