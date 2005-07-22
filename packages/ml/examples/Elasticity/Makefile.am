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

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_example_elasticity \
                  ml_read_elas

ml_example_elasticity_SOURCES      = $(srcdir)/ml_example_elasticity.c
ml_example_elasticity_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_example_elasticity_CFLAGS       = $(EXAMPLE_INCLUDES)
ml_example_elasticity_LDADD        = $(EXAMPLE_LIBS)
ml_example_elasticity_LINK         = $(ML_LINK)

ml_read_elas_SOURCES      = $(srcdir)/ml_read_elas.c
ml_read_elas_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_read_elas_CFLAGS       = $(EXAMPLE_INCLUDES)
ml_read_elas_LDADD        = $(EXAMPLE_LIBS)
ml_read_elas_LINK         = $(ML_LINK)

