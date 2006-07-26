#@HEADER

# include definitions of ML_INCLUDES and ML_LIBS
include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_AztecMsr   \
                  ml_EpetraVbr \
                  ml_MatrixFree

ml_AztecMsr_SOURCES      = $(srcdir)/ml_AztecMsr.cpp
ml_AztecMsr_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_AztecMsr_CXXFLAGS     = ${EXAMPLE_INCLUDES} 
ml_AztecMsr_LDADD        = ${EXAMPLE_LIBS} 

ml_EpetraVbr_SOURCES      = $(srcdir)/ml_EpetraVbr.cpp
ml_EpetraVbr_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_EpetraVbr_CXXFLAGS     = ${EXAMPLE_INCLUDES}
ml_EpetraVbr_LDADD        = ${EXAMPLE_LIBS}

ml_MatrixFree_SOURCES      = $(srcdir)/ml_MatrixFree.cpp
ml_MatrixFree_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_MatrixFree_CXXFLAGS     = ${EXAMPLE_INCLUDES}
ml_MatrixFree_LDADD        = ${EXAMPLE_LIBS}
