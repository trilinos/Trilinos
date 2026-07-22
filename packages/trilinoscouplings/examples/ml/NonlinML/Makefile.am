
# include definitions of ML_INCLUDES and ML_LIBS
include $(top_builddir)/Makefile.export.trilinoscouplings

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_nox_1Delasticity_example

ml_nox_1Delasticity_example_SOURCES = $(srcdir)/ml_nox_1Delasticity_example.cpp \
                                      $(srcdir)/Basis.cpp $(srcdir)/Basis.H \
                                      $(srcdir)/FiniteElementProblem.cpp $(srcdir)/FiniteElementProblem.H \
                                      $(srcdir)/Problem_Interface.cpp $(srcdir)/Problem_Interface.H
ml_nox_1Delasticity_example_CXXFLAGS = $(NOX_INCLUDES) $(TRILINOSCOUPLINGS_INCLUDES)
ml_nox_1Delasticity_example_LDADD = $(NOX_LIBS) $(TRILINOSCOUPLINGS_LIBS)
#ml_nox_1Delasticity_example_LINK = ${mllink}
ml_nox_1Delasticity_example_DEPENDENCIES = $(NOX_INCLUDES) $(TRILINOSCOUPLINGS_LIBS)
