#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
ml_common_include_locations = $(ML_EXPORT_INCLUDES)
ml_libs = $(ML_EXPORT_LIBS) $(LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
mllink = $(CXX) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_viz
ml_viz_SOURCES = $(srcdir)/ml_viz.cpp
ml_viz_CPPFLAGS = ${ml_common_include_locations}
ml_viz_CXXFLAGS = ${ml_common_include_locations}
ml_viz_LDADD = ${ml_libs} 
ml_viz_LINK = ${mllink}
ml_viz_DEPENDENCIES = $(top_builddir)/src/libml.a
