#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_example_elasticity \
                  ml_read_elas

ml_example_elasticity_SOURCES = $(srcdir)/ml_example_elasticity.c
ml_example_elasticity_CFLAGS = ${ml_common_include_locations}
ml_example_elasticity_LDADD = ${ml_libs}
ml_example_elasticity_LINK = ${mllink}
ml_example_elasticity_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_read_elas_SOURCES = $(srcdir)/ml_read_elas.c
ml_read_elas_CFLAGS = ${ml_common_include_locations}
ml_read_elas_LDADD = ${ml_libs}
ml_read_elas_LINK = ${mllink}
ml_read_elas_DEPENDENCIES = $(top_builddir)/src/libml.a
