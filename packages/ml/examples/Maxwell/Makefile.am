#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_maxwell \
                  ml_simple_max

ml_maxwell_SOURCES = $(srcdir)/ml_maxwell.cpp
ml_maxwell_CPPFLAGS = ${ml_common_include_locations}
ml_maxwell_CXXFLAGS = ${ml_common_include_locations}
ml_maxwell_LDADD = ${ml_libs}
ml_maxwell_LINK = ${mllink}
ml_maxwell_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_simple_max_SOURCES = $(srcdir)/ml_simple_max.c
ml_simple_max_CFLAGS = ${ml_common_include_locations}
ml_simple_max_LDADD = ${ml_libs}
ml_simple_max_LINK = ${mllink}
ml_simple_max_DEPENDENCIES = $(top_builddir)/src/libml.a
