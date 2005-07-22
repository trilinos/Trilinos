#@HEADER

include $(top_builddir)/Makefile.export.ml

if USING_GNUMAKE
EXAMPLE_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ML_LIBS) $(ML_EXTRA_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ML_INCLUDES) $(ML_EXTRA_INCLUDES))
else
EXAMPLE_LIBS     = $(ML_LIBS) $(ML_EXTRA_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES) $(ML_EXTRA_INCLUDES)
endif

AM_CXXFLAGS = $(EXAMPLE_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_aztec_simple          \
                  ml_operator              \
                  ml_preconditioner        \
		  ml_readex

ml_aztec_simple_SOURCES      = $(srcdir)/ml_aztec_simple.c
ml_aztec_simple_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_aztec_simple_CFLAGS       = $(EXAMPLE_INCLUDES)
ml_aztec_simple_LDADD        = $(EXAMPLE_LIBS)
ml_aztec_simple_LINK         = $(ML_LINK)

ml_operator_SOURCES      = $(srcdir)/ml_operator.cpp
ml_operator_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_operator_LDADD        = $(EXAMPLE_LIBS)

ml_preconditioner_SOURCES      = $(srcdir)/ml_preconditioner.cpp
ml_preconditioner_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_preconditioner_LDADD        = $(EXAMPLE_LIBS)

ml_readex_SOURCES      = $(srcdir)/ml_readex.c
ml_readex_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_readex_CFLAGS       = $(EXAMPLE_INCLUDES)
ml_readex_LDADD        = $(EXAMPLE_LIBS)
ml_readex_LINK         = $(ML_LINK)
