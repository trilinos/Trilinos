#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_read_matrix_AIJ
                  

ml_read_matrix_AIJ_SOURCES= $(srcdir)/ml_read_matrix_AIJ.cpp
ml_read_matrix_AIJ_CPPFLAGS = ${ml_common_include_locations}
ml_read_matrix_AIJ_CXXFLAGS = ${ml_common_include_locations}
ml_read_matrix_AIJ_LDADD = ${ml_libs}
ml_read_matrix_AIJ_LINK = ${mllink}
ml_read_matrix_AIJ_DEPENDENCIES = $(top_builddir)/src/libml.a
