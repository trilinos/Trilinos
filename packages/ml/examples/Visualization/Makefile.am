#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_viz
ml_viz_SOURCES = $(srcdir)/ml_viz.cpp
ml_viz_CPPFLAGS = ${ml_common_include_locations}
ml_viz_CXXFLAGS = ${ml_common_include_locations}
ml_viz_LDADD = ${ml_libs} 
ml_viz_LINK = ${mllink}
ml_viz_DEPENDENCIES = $(top_builddir)/src/libml.a
