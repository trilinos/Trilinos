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
                                       -I$(top_srcdir)/../nox/src \
                                       -I$(top_srcdir)/../nox/src-epetra \
                                       -I$(top_srcdir)/../epetraext/src/coloring \
                                       -I$(top_srcdir)/../epetraext/src/transform \
                                       -I$(top_builddir)/../nox/src
ml_nox_1Delasticity_example_LDADD = ${ml_libs} \
                                    $(top_builddir)/src/libml.a \
                                    $(top_builddir)/../nox/src-epetra/libnoxepetra.a \
                                    $(top_builddir)/../nox/src/libnox.a \
                                    $(top_builddir)/../aztecoo/src/libaztecoo.a \
                                    $(top_builddir)/../ifpack/src/libifpack.a \
                                    $(top_builddir)/../amesos/src/libamesos.a \
                                    $(top_builddir)/../teuchos/src/libteuchos.a \
                                    $(top_builddir)/../epetraext/src/libepetraext.a \
                                    $(top_builddir)/../epetra/src/libepetra.a \
                                    $(top_builddir)/src/libml.a

ml_nox_1Delasticity_example_LINK = ${mllink}
ml_nox_1Delasticity_example_DEPENDENCIES = $(top_builddir)/src/libml.a \
                                           $(top_builddir)/../nox/src-epetra/libnoxepetra.a \
                                           $(top_builddir)/../nox/src/libnox.a \
                                           $(top_builddir)/../aztecoo/src/libaztecoo.a \
                                           $(top_builddir)/../ifpack/src/libifpack.a \
                                           $(top_builddir)/../amesos/src/libamesos.a \
                                           $(top_builddir)/../teuchos/src/libteuchos.a \
                                           $(top_builddir)/../epetraext/src/libepetraext.a \
                                           $(top_builddir)/../epetra/src/libepetra.a \
                                           $(top_builddir)/src/libml.a
