#@HEADER

include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_example1d  \
                  ml_example1dGS    \
                  mlguide   \
                  mlguide_par

CCLD = $(CXX)

mlguide_SOURCES      = $(srcdir)/mlguide.c
mlguide_DEPENDENCIES = $(top_builddir)/src/libml.a
mlguide_LDADD        = ${EXAMPLE_LIBS}
mlguide_CFLAGS       = ${EXAMPLE_INCLUDES}

mlguide_par_SOURCES      = $(srcdir)/mlguide_par.c
mlguide_par_DEPENDENCIES = $(top_builddir)/src/libml.a
mlguide_par_LDADD        = ${EXAMPLE_LIBS}
mlguide_par_CFLAGS       = ${EXAMPLE_INCLUDES}

ml_example1d_SOURCES      = $(srcdir)/ml_example1d.c
ml_example1d_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_example1d_LDADD        = ${EXAMPLE_LIBS}
ml_example1d_CFLAGS       = ${EXAMPLE_INCLUDES}

ml_example1dGS_SOURCES      = $(srcdir)/ml_example1dGS.c
ml_example1dGS_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_example1dGS_LDADD        = ${EXAMPLE_LIBS}
ml_example1dGS_CFLAGS       = ${EXAMPLE_INCLUDES}
