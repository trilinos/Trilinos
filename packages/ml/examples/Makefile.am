#@HEADER
#@HEADER

EXEEXT = .exe

#EXAMPLE = ml_read_maxwell

# List the names of all exe's in this directory here
#noinst_PROGRAMS = ml_read_maxwell
noinst_PROGRAMS = mlguide

# One SOURCES primary is needed for EACH exe listed above.  (Notice for all 
# of the  following primaries that the pattern is 
# "(Exe_name)_(PRIMARY_NAME) = ").
# List all source files needed to compile the exe.
ml_read_maxwell_SOURCES = $(srcdir)/ml_read_maxwell.cpp
mlguide_SOURCES = $(srcdir)/mlguide.c

# One DEPENDENCIES primary is needed for each exe listed above.  List any
# libraries that the exe depends on.  (This primary is not absolutely
# essential.)
ml_read_maxwell_DEPENDENCIES = $(top_builddir)/../aztecoo/src/libaztecoo.a \
                       $(top_builddir)/../epetra/src/libepetra.a
mlguide_DEPENDENCIES = 
                       

ml_common_include_locations = -I$(top_srcdir)/src/Include          \
                   -I$(top_srcdir)/src/Coarsen                  \
                   -I$(top_srcdir)/src/Comm                  \
                   -I$(top_srcdir)/src/FEGrid                  \
                   -I$(top_srcdir)/src/Include                  \
                   -I$(top_srcdir)/src/Krylov                  \
                   -I$(top_srcdir)/src/Main                  \
                   -I$(top_srcdir)/src/Operator                  \
                   -I$(top_srcdir)/src/Smoother                  \
                   -I$(top_srcdir)/src/Utils                  \
                   -I$(top_srcdir)/../aztecoo/src            \
                   -I$(top_builddir)/../aztecoo/src       

# One CXXFLAGS primary can be defined for each exe listed above.  Include 
# any -I's necessary for building the exe.
ml_read_maxwell_CXXFLAGS = ${ml_common_include_locations}
mlguide_CFLAGS = ${ml_common_include_locations}

# One LDADD primary can be defined for each exe listed above.  Include any
# LDFLAGS  necessary for building the exe.
ml_read_maxwell_LDADD = $(top_builddir)/../aztecoo/src/libaztecoo.a        \
                $(top_builddir)/src/libml.a        \
                $(top_builddir)/../epetra/src/libepetra.a         \
                $(LAPACK_LIBS)                            \
                $(BLAS_LIBS)                            \
                $(LIBS) $(FLIBS)
#
mlguide_LDADD = $(top_builddir)/src/libml.a        \
                $(LAPACK_LIBS)                            \
                $(BLAS_LIBS)                            \
                $(LIBS) $(FLIBS)
