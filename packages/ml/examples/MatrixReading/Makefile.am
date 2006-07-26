#@HEADER

include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_read_matrix_AIJ

ml_read_matrix_AIJ_SOURCES      = $(srcdir)/ml_read_matrix_AIJ.cpp
ml_read_matrix_AIJ_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_read_matrix_AIJ_CXXFLAGS     = ${EXAMPLE_INCLUDES}
ml_read_matrix_AIJ_LDADD        = ${EXAMPLE_LIBS}
