#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_example1d  \
                  ml_example1dGS    \
                  mlguide   \
                  mlguide_par


###############
# sources
###############

mlguide_SOURCES = $(srcdir)/mlguide.c
mlguide_CFLAGS = ${ml_common_include_locations}
mlguide_LDADD = ${ml_libs}
mlguide_LINK = ${mllink}
mlguide_DEPENDENCIES = $(top_builddir)/src/libml.a

mlguide_par_SOURCES = $(srcdir)/mlguide_par.c
mlguide_par_CFLAGS = ${ml_common_include_locations}
mlguide_par_LDADD = ${ml_libs}
mlguide_par_LINK = ${mllink}
mlguide_par_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_example1d_SOURCES = $(srcdir)/ml_example1d.c
ml_example1d_CFLAGS = ${ml_common_include_locations}
ml_example1d_LDADD = ${ml_libs}
ml_example1d_LINK = ${mllink}
ml_example1d_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_example1dGS_SOURCES = $(srcdir)/ml_example1dGS.c
ml_example1dGS_CFLAGS = ${ml_common_include_locations}
ml_example1dGS_LDADD = ${ml_libs}
ml_example1dGS_LINK = ${mllink}
ml_example1dGS_DEPENDENCIES = $(top_builddir)/src/libml.a
