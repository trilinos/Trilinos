#@HEADER

include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_viz
ml_viz_SOURCES      = $(srcdir)/ml_viz.cpp
ml_viz_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_viz_CXXFLAGS     = ${EXAMPLE_INCLUDES} 
ml_viz_LDADD        = ${EXAMPLE_LIBS} 
