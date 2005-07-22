#@HEADER

include $(top_builddir)/Makefile.export.ml

if USING_GNUMAKE
EXAMPLE_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ML_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ML_INCLUDES))
else
EXAMPLE_LIBS     = $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES)
endif

AM_CPPFLAGS = $(EXAMPLE_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_read_matrix_AIJ

ml_read_matrix_AIJ_SOURCES      = $(srcdir)/ml_read_matrix_AIJ.cpp
ml_read_matrix_AIJ_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_read_matrix_AIJ_LDADD        = ${EXAMPLE_LIBS}
