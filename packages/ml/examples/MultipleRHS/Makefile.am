#@HEADER

include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_LIBS) 
EXAMPLE_INCLUDES = $(ML_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_multiple_rhs

ml_multiple_rhs_SOURCES      = $(srcdir)/ml_multiple_rhs.cpp
ml_multiple_rhs_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_multiple_rhs_CXXFLAGS     = $(EXAMPLE_INCLUDES)
ml_multiple_rhs_LDADD        = $(EXAMPLE_LIBS)
