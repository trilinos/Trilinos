#@HEADER

include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_aztec_simple   \
                  ml_operator       \
                  ml_preconditioner \
                  ml_user_smoothing \
                  ml_readex \
                  GetNecessaryFiles

CCLD = $(CXX)

ml_aztec_simple_SOURCES      = $(srcdir)/ml_aztec_simple.c
ml_aztec_simple_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_aztec_simple_CFLAGS       = $(EXAMPLE_INCLUDES)
ml_aztec_simple_LDADD        = $(EXAMPLE_LIBS)

ml_operator_SOURCES      = $(srcdir)/ml_operator.cpp
ml_operator_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_operator_CXXFLAGS     = $(EXAMPLE_INCLUDES)
ml_operator_LDADD        = $(EXAMPLE_LIBS)

ml_preconditioner_SOURCES      = $(srcdir)/ml_preconditioner.cpp
ml_preconditioner_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_preconditioner_CXXFLAGS     = $(EXAMPLE_INCLUDES)
ml_preconditioner_LDADD        = $(EXAMPLE_LIBS)

ml_readex_SOURCES      = $(srcdir)/ml_readex.c
ml_readex_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_readex_CFLAGS       = $(EXAMPLE_INCLUDES)
ml_readex_LDADD        = $(EXAMPLE_LIBS)

ml_user_smoothing_SOURCES      = $(srcdir)/ml_user_smoothing.cpp
ml_user_smoothing_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_user_smoothing_CXXFLAGS     = $(EXAMPLE_INCLUDES)
ml_user_smoothing_LDADD        = $(EXAMPLE_LIBS)

GetNecessaryFiles_SOURCES= $(top_srcdir)/examples/BasicExamples
#copy these files over to the build area
GetNecessaryFiles$(EXEEXT): $(GetNecessaryFiles_SOURCES)
	@if test ! -e ./.data; then cp -f $(srcdir)/.data.gz .; gunzip -f ./.data.gz ; fi
