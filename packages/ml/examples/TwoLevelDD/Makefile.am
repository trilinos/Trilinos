#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_2level_DD

ml_2level_DD_SOURCES = $(srcdir)/ml_2level_DD.cpp
ml_2level_DD_CPPFLAGS = ${ml_common_include_locations}
ml_2level_DD_CXXFLAGS = ${ml_common_include_locations}
ml_2level_DD_LDADD = ${ml_libs} 
ml_2level_DD_LINK = ${mllink}
ml_2level_DD_DEPENDENCIES = $(top_builddir)/src/libml.a
