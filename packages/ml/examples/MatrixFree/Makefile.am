#@HEADER

# include definitions of ML_INCLUDES and ML_LIBS
# ML_EXTRA_LIBS and ML_EXTRA_INCLUDES contain additional libraries and
# header files required to compile and link some of the ML examples.
include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_1D

ml_1D_SOURCES      = $(srcdir)/ml_1D.cpp
ml_1D_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_1D_CXXFLAGS     = $(EXAMPLE_INCLUDES)
ml_1D_LDADD        = $(EXAMPLE_LIBS)
