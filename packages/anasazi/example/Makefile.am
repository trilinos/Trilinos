#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = BlockArnoldiPetraEx

BlockArnoldiPetraEx_SOURCES = $(srcdir)/BlockArnoldiPetraEx.cpp

#np# Specify locations of include files within Trilinos, but outside np
#np# Add packages that np is dependent on; remove epetra if np does not need it

AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src

#np# - List libraries that np is dependent on.  The DEPENDENCIES primary is
#np# optional.
BlockArnoldiPetraEx_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/utils/libanasaziutils.a 

BlockArnoldiPetraEx_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/utils \
				-I$(top_srcdir)/../epetra/src

#np# List libraries that np is dependent on.
BlockArnoldiPetraEx_LDADD = $(top_builddir)/utils/libanasaziutils.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose examples are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# example1/ rather than this file.
# SUBDIRS = example1
SUBDIRS = 
