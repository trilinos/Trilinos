#@HEADER

include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_example_elasticity \
                  ml_read_elas

CCLD = $(CXX)

ml_example_elasticity_SOURCES      = $(srcdir)/ml_example_elasticity.c
ml_example_elasticity_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_example_elasticity_CFLAGS       = $(EXAMPLE_INCLUDES)
ml_example_elasticity_LDADD        = $(EXAMPLE_LIBS)

ml_read_elas_SOURCES      = $(srcdir)/ml_read_elas.c
ml_read_elas_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_read_elas_CFLAGS       = $(EXAMPLE_INCLUDES)
ml_read_elas_LDADD        = $(EXAMPLE_LIBS)

