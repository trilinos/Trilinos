#@HEADER

include $(top_builddir)/Makefile.export.ml

if USING_GNUMAKE
EXAMPLE_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ML_LIBS) $(ML_EXTRA_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ML_INCLUDES) $(ML_EXTRA_INCLUDES))
else
EXAMPLE_LIBS     = $(ML_LIBS)     $(ML_EXTRA_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES) $(ML_EXTRA_INCLUDES)
endif

AM_CPPFLAGS = $(EXAMPLE_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_2level_DD

ml_2level_DD_SOURCES      = $(srcdir)/ml_2level_DD.cpp
ml_2level_DD_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_2level_DD_LDADD        = ${EXAMPLE_LIBS} 
