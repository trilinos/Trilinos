#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
ml_common_include_locations = $(ML_EXPORT_INCLUDES)
ml_libs = $(ML_EXPORT_LIBS) $(LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
mllink = $(CXX) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_user_aggr \
                  ml_filtering \
                  ml_analyze 

ml_user_aggr_SOURCES = $(srcdir)/ml_user_aggr.cpp
ml_user_aggr_CPPFLAGS = ${ml_common_include_locations}
ml_user_aggr_CXXFLAGS = ${ml_common_include_locations}
ml_user_aggr_LDADD = ${ml_libs} 
ml_user_aggr_LINK = ${mllink}
ml_user_aggr_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_filtering_SOURCES = $(srcdir)/ml_filtering.cpp
ml_filtering_CPPFLAGS = ${ml_common_include_locations}
ml_filtering_CXXFLAGS = ${ml_common_include_locations}
ml_filtering_LDADD = ${ml_libs} 
ml_filtering_LINK = ${mllink}
ml_filtering_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_analyze_SOURCES = $(srcdir)/ml_analyze.cpp
ml_analyze_CPPFLAGS = ${ml_common_include_locations}
ml_analyze_CXXFLAGS = ${ml_common_include_locations}
ml_analyze_LDADD = ${ml_libs} 
ml_analyze_LINK = ${mllink}
ml_analyze_DEPENDENCIES = $(top_builddir)/src/libml.a
