#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_example2d \
                  ml_example3d    \
                  ml_read_elas    \
                  ml_readfile    \
                  ml_read_salsa    \
                  ml_recirc    

ml_example2d_SOURCES = $(srcdir)/ml_example2d.c
ml_example2d_CFLAGS = ${ml_common_include_locations}
ml_example2d_LDADD = ${ml_libs}
ml_example2d_LINK = ${mllink}
ml_example2d_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_example3d_SOURCES = $(srcdir)/ml_example3d.c
ml_example3d_CFLAGS = ${ml_common_include_locations}
ml_example3d_LDADD = ${ml_libs}
ml_example3d_LINK = ${mllink}
ml_example3d_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_read_elas_SOURCES = $(srcdir)/ml_read_elas.c
ml_read_elas_CFLAGS = ${ml_common_include_locations}
ml_read_elas_LDADD = ${ml_libs}
ml_read_elas_LINK = ${mllink}
ml_read_elas_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_readfile_SOURCES = $(srcdir)/ml_readfile.c
ml_readfile_CFLAGS = ${ml_common_include_locations}
ml_readfile_LDADD = ${ml_libs}
ml_readfile_LINK = ${mllink}
ml_readfile_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_read_salsa_SOURCES = $(srcdir)/ml_read_salsa.c
ml_read_salsa_CFLAGS = ${ml_common_include_locations}
ml_read_salsa_LDADD = ${ml_libs}
ml_read_salsa_LINK = ${mllink}
ml_read_salsa_DEPENDENCIES = $(top_builddir)/src/libml.a

ml_recirc_SOURCES = $(srcdir)/ml_recirc.c
ml_recirc_CFLAGS = ${ml_common_include_locations}
ml_recirc_LDADD = ${ml_libs}
ml_recirc_LINK = ${mllink}
ml_recirc_DEPENDENCIES = $(top_builddir)/src/libml.a
