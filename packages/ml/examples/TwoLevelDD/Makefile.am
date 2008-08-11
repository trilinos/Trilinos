#@HEADER

include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_2level_DD

ml_2level_DD_SOURCES      = $(srcdir)/ml_2level_DD.cpp
ml_2level_DD_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_2level_DD_CXXFLAGS     = ${EXAMPLE_INCLUDES} 
ml_2level_DD_LDADD        = ${EXAMPLE_LIBS} 
