#@HEADER

include $(top_builddir)/Makefile.export.ml

if USING_GNUMAKE
EXAMPLE_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ML_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ML_INCLUDES))
else
EXAMPLE_LIBS     = $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES)
endif

EXEEXT = .exe

AM_CXXFLAGS = ${EXAMPLE_INCLUDES}

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_maxwell \
                  ml_simple_max

ml_maxwell_SOURCES      = $(srcdir)/ml_maxwell.cpp
ml_maxwell_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_maxwell_LDADD        = ${EXAMPLE_LIBS}

ml_simple_max_SOURCES      = $(srcdir)/ml_simple_max.c
ml_simple_max_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_simple_max_CFLAGS       = ${EXAMPLE_INCLUDES}
ml_simple_max_LDADD        = ${EXAMPLE_LIBS}
ml_simple_max_LINK         = ${ML_LINK}
