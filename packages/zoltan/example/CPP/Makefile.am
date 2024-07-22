
# The following line helps the test harness recover from build errors.

all-local:

include $(top_builddir)/Makefile.export.zoltan

EXEEXT = .exe

noinst_PROGRAMS = exampleBLOCK

exampleBLOCK_SOURCES      = exampleBLOCK.cpp
exampleBLOCK_DEPENDENCIES = $(ZOLTAN_DEPS)
exampleBLOCK_CXXFLAGS     = $(ZOLTAN_INCLUDES)
exampleBLOCK_LDADD        = $(ZOLTAN_LIBS)

EXTRA_DIST = classicMakefile objects.txt
