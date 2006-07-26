#@HEADER
#@HEADER

# include definitions of ML_INCLUDES and ML_LIBS
include $(top_builddir)/Makefile.export.ml

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_nox_1Delasticity_example

ml_nox_1Delasticity_example_SOURCES = $(srcdir)/ml_nox_1Delasticity_example.cpp \
                                      $(srcdir)/Basis.cpp $(srcdir)/Basis.H \
                                      $(srcdir)/FiniteElementProblem.cpp $(srcdir)/FiniteElementProblem.H \
                                      $(srcdir)/Problem_Interface.cpp $(srcdir)/Problem_Interface.H
ml_nox_1Delasticity_example_CXXFLAGS = -I$(top_srcdir)/src/NonlinML \
                                       $(ML_INCLUDES)
ml_nox_1Delasticity_example_LDADD = -L$(ML_BUILD_DIRECTORY)/src -lml \
                                    $(NOX_LIBS) $(ML_LIBS)
#ml_nox_1Delasticity_example_LINK = ${mllink}
ml_nox_1Delasticity_example_DEPENDENCIES = $(top_builddir)/src/libml.a \
                                           $(ML_LIBS) \
                                           $(top_builddir)/src/libml.a
