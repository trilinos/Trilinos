#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_AztecMsr   \
                  ml_EpetraVbr \
                  ml_MatrixFree

ml_AztecMsr_SOURCES = $(srcdir)/ml_AztecMsr.cpp
ml_AztecMsr_CPPFLAGS = ${ml_common_include_locations}
ml_AztecMsr_CXXFLAGS = ${ml_common_include_locations}
ml_AztecMsr_LDADD = ${ml_libs} 
ml_AztecMsr_LINK = ${mllink}
ml_AztecMsr_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_EpetraVbr_SOURCES= $(srcdir)/ml_EpetraVbr.cpp
ml_EpetraVbr_CPPFLAGS = ${ml_common_include_locations}
ml_EpetraVbr_CXXFLAGS = ${ml_common_include_locations}
ml_EpetraVbr_LDADD = ${ml_libs}
ml_EpetraVbr_LINK = ${mllink}
ml_EpetraVbr_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_MatrixFree_SOURCES= $(srcdir)/ml_MatrixFree.cpp
ml_MatrixFree_CPPFLAGS = ${ml_common_include_locations}
ml_MatrixFree_CXXFLAGS = ${ml_common_include_locations}
ml_MatrixFree_LDADD = ${ml_libs}
ml_MatrixFree_LINK = ${mllink}
ml_MatrixFree_DEPENDENCIES = $(top_builddir)/src/libml.a
