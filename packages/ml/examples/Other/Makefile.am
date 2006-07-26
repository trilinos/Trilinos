#@HEADER

include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_LIBS)     $(ML_EXTRA_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES) $(ML_EXTRA_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_example2d   \
                  ml_example3d   \
                  ml_read_elas   \
                  ml_readfile    \
                  ml_read_salsa  \
                  ml_recirc    

CCLD = $(CXX)

ml_example2d_SOURCES      = $(srcdir)/ml_example2d.c
ml_example2d_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_example2d_LDADD        = ${EXAMPLE_LIBS}
ml_example2d_CFLAGS       = ${EXAMPLE_INCLUDES}

ml_example3d_SOURCES      = $(srcdir)/ml_example3d.c
ml_example3d_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_example3d_LDADD        = ${EXAMPLE_LIBS}
ml_example3d_CFLAGS       = ${EXAMPLE_INCLUDES}

ml_read_elas_SOURCES      = $(srcdir)/ml_read_elas.c
ml_read_elas_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_read_elas_LDADD        = ${EXAMPLE_LIBS}
ml_read_elas_CFLAGS       = ${EXAMPLE_INCLUDES}

ml_readfile_SOURCES      = $(srcdir)/ml_readfile.c
ml_readfile_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_readfile_LDADD        = ${EXAMPLE_LIBS}
ml_readfile_CFLAGS       = ${EXAMPLE_INCLUDES}

ml_read_salsa_SOURCES      = $(srcdir)/ml_read_salsa.c
ml_read_salsa_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_read_salsa_LDADD        = ${EXAMPLE_LIBS}
ml_read_salsa_CFLAGS       = ${EXAMPLE_INCLUDES}

ml_recirc_SOURCES      = $(srcdir)/ml_recirc.c
ml_recirc_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_recirc_LDADD        = ${EXAMPLE_LIBS}
ml_recirc_CFLAGS       = ${EXAMPLE_INCLUDES}
