
include $(top_builddir)/Makefile.export.trilinoscouplings

TEST_LIBS     = $(TRILINOSCOUPLINGS_LIBS)
TEST_INCLUDES = $(TRILINOSCOUPLINGS_INCLUDES)

EXEEXT = .exe

if BUILD_IFPACK_ML

noinst_PROGRAMS = DDwithMLSubsolves

DDwithMLSubsolves_SOURCES      = $(srcdir)/cxx_main.cpp
DDwithMLSubsolves_DEPENDENCIES = $(top_builddir)/src/libml.a
DDwithMLSubsolves_CXXFLAGS     = ${TEST_INCLUDES}
DDwithMLSubsolves_LDADD        = ${TEST_LIBS}

endif

