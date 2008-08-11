#@HEADER

include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_maxwell \
                  ml_simple_max \
                  ml_read_maxwell

CCLD = $(CXX)

ml_maxwell_SOURCES      = $(srcdir)/ml_maxwell.cpp
ml_maxwell_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_maxwell_CXXFLAGS     = ${EXAMPLE_INCLUDES}
ml_maxwell_LDADD        = ${EXAMPLE_LIBS}

ml_simple_max_SOURCES      = $(srcdir)/ml_simple_max.c
ml_simple_max_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_simple_max_CFLAGS       = ${EXAMPLE_INCLUDES}
ml_simple_max_LDADD        = ${EXAMPLE_LIBS}

ml_read_maxwell_SOURCES      = $(srcdir)/ml_read_maxwell.cpp
ml_read_maxwell_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_read_maxwell_CXXFLAGS     = ${EXAMPLE_INCLUDES}
ml_read_maxwell_LDADD        = ${EXAMPLE_LIBS}
