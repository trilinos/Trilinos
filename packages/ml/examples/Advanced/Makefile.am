#@HEADER

# include definitions of ML_INCLUDES and ML_LIBS
# ML_EXTRA_LIBS and ML_EXTRA_INCLUDES contain additional libraries and
# header files required to compile and link some of the ML examples.
include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_user_aggr \
                  ml_filtering \
                  ml_analyze   \
		  ml_with_ifpack

ml_user_aggr_SOURCES      = $(srcdir)/ml_user_aggr.cpp
ml_user_aggr_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_user_aggr_CXXFLAGS     = $(EXAMPLE_INCLUDES)
ml_user_aggr_LDADD        = $(EXAMPLE_LIBS)

ml_filtering_SOURCES      = $(srcdir)/ml_filtering.cpp
ml_filtering_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_filtering_CXXFLAGS     = $(EXAMPLE_INCLUDES)
ml_filtering_LDADD        = $(EXAMPLE_LIBS)

ml_analyze_SOURCES      = $(srcdir)/ml_analyze.cpp
ml_analyze_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_analyze_CXXFLAGS     = $(EXAMPLE_INCLUDES)
ml_analyze_LDADD        = $(EXAMPLE_LIBS)

ml_with_ifpack_SOURCES      = $(srcdir)/ml_with_ifpack.cpp
ml_with_ifpack_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_with_ifpack_CXXFLAGS     = $(EXAMPLE_INCLUDES)
ml_with_ifpack_LDADD        = $(EXAMPLE_LIBS)
