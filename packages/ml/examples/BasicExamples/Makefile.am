#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
ml_common_include_locations = $(ML_EXPORT_INCLUDES)
ml_libs = $(ML_EXPORT_LIBS) $(LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
mllink = $(CXX) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_aztec_simple          \
                  ml_operator              \
                  ml_preconditioner        \
                  ml_readex

ml_aztec_simple_SOURCES = $(srcdir)/ml_aztec_simple.c
ml_aztec_simple_CFLAGS = ${ml_common_include_locations}
ml_aztec_simple_LDADD = ${ml_libs} 
ml_aztec_simple_LINK = ${mllink}
ml_aztec_simple_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_operator_SOURCES = $(srcdir)/ml_operator.cpp
ml_operator_CPPFLAGS = ${ml_common_include_locations}
ml_operator_CXXFLAGS = ${ml_common_include_locations}
ml_operator_LDADD = ${ml_libs} 
ml_operator_LINK = ${mllink}
ml_operator_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_preconditioner_SOURCES = $(srcdir)/ml_preconditioner.cpp
ml_preconditioner_CPPFLAGS = ${ml_common_include_locations}
ml_preconditioner_CXXFLAGS = ${ml_common_include_locations}
ml_preconditioner_LDADD = ${ml_libs} 
ml_preconditioner_LINK = ${mllink}
ml_preconditioner_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_readex_SOURCES = $(srcdir)/ml_readex.c
ml_readex_CFLAGS = ${ml_common_include_locations}
ml_readex_LDADD = ${ml_libs} 
ml_readex_LINK = ${mllink}
ml_readex_DEPENDENCIES = $(top_builddir)/src/libml.a

