#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_multiple_rhs

ml_multiple_rhs_SOURCES = $(srcdir)/ml_multiple_rhs.cpp
ml_multiple_rhs_CPPFLAGS = ${ml_common_include_locations}
ml_multiple_rhs_CXXFLAGS = ${ml_common_include_locations}
ml_multiple_rhs_LDADD = ${ml_libs}
ml_multiple_rhs_LINK = ${mllink}
ml_multiple_rhs_DEPENDENCIES = $(top_builddir)/src/libml.a
