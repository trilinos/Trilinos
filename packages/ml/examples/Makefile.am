#@HEADER
#@HEADER
#np# This file is an example of a typical Makefile.am to build an example exe. 

EXEEXT = .exe

EXAMPLE = ml_read_maxwell

#np# List the names of all exe's in this directory here
noinst_PROGRAMS = ml_read_maxwell

#np# One SOURCES primary is needed for EACH exe listed above.  (Notice for all 
#np# of the  following primaries that the pattern is 
#np# "(Exe_name)_(PRIMARY_NAME) = ").
#np# List all source files needed to compile the exe.
ml_read_maxwell_SOURCES = $(srcdir)/ml_read_maxwell.c

#np# One DEPENDENCIES primary is needed for each exe listed above.  List any
#np# libraries that the exe depends on.  (This primary is not absolutely
#np# essential.)
ml_read_maxwell_DEPENDENCIES = $(top_builddir)/../aztecoo/src/libaztecoo.a \
                       $(top_builddir)/../epetra/src/libepetra.a

#np# One CXXFLAGS primary can be defined for each exe listed above.  Include 
#np# any -I's necessary for building the exe.
ml_read_maxwell_CFLAGS = -I$(top_srcdir)/src/Include          \
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

#np# One LDADD primary can be defined for each exe listed above.  Include any
#np# LDFLAGS  necessary for building the exe.
ml_read_maxwell_LDADD = $(top_builddir)/../aztecoo/src/libaztecoo.a        \
                $(top_builddir)/src/libml.a        \
                $(top_builddir)/../epetra/src/libepetra.a         \
                $(LAPACK_LIBS)                            \
                $(BLAS_LIBS)                            \
                $(LIBS) $(FLIBS)
