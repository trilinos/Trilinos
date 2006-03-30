#@HEADER
#@HEADER

# include definitions of ML_INCLUDES and ML_LIBS
include $(top_builddir)/Makefile.export

EXAMPLE_LIBS     = $(ML_LIBS)     $(ML_EXTRA_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES) $(ML_EXTRA_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_nox_1Delasticity_example

ml_nox_1Delasticity_example_SOURCES      = $(srcdir)/ml_nox_1Delasticity_example.cpp \
                                           $(srcdir)/Basis.cpp $(srcdir)/Basis.H \
                                           $(srcdir)/FiniteElementProblem.cpp $(srcdir)/FiniteElementProblem.H \
                                           $(srcdir)/Problem_Interface.cpp $(srcdir)/Problem_Interface.H
ml_nox_1Delasticity_example_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_nox_1Delasticity_example_CXXFLAGS     = $(EXAMPLE_INCLUDES)
ml_nox_1Delasticity_example_LDADD        = ${EXAMPLE_LIBS} $(NOX_LIBS) 
