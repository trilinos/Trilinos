#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_nox_1Delasticity_example

ml_nox_1Delasticity_example_SOURCES = $(srcdir)/ml_nox_1Delasticity_example.cpp \
                                      $(srcdir)/Basis.cpp $(srcdir)/Basis.H \
                                      $(srcdir)/FiniteElementProblem.cpp $(srcdir)/FiniteElementProblem.H \
                                      $(srcdir)/Problem_Interface.cpp $(srcdir)/Problem_Interface.H
ml_nox_1Delasticity_example_CPPFLAGS = ${ml_common_include_locations}
ml_nox_1Delasticity_example_CXXFLAGS = ${ml_common_include_locations} \
                                       -I$(top_srcdir)/src/NonlinML \
                                       -I$(ML_EXPORT_INCLUDES)
ml_nox_1Delasticity_example_LDADD = ${ml_libs} \
                                    $(top_builddir)/src/libml.a \
                                    $(ML_EXPORT_LIBS) \
                                    $(top_builddir)/src/libml.a
ml_nox_1Delasticity_example_LINK = ${mllink}
ml_nox_1Delasticity_example_DEPENDENCIES = $(top_builddir)/src/libml.a \
                                           $(ML_EXPORT_LIBS) \
                                           $(top_builddir)/src/libml.a
