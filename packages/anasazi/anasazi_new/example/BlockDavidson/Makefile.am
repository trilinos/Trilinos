#@HEADER
#@HEADER

EXEEXT = .exe

if BUILD_BELOS
if BUILD_IFPACK
ANASAZI_BIRAM_EX_SI_IFPACK = BlockDavidsonEpetraExGen
endif
endif

noinst_PROGRAMS = BlockDavidsonEpetraEx $(ANASAZI_BIRAM_EX_HB) $(ANASAZI_BIRAM_EX_SI_IFPACK) 

BlockDavidsonEpetraEx_SOURCES = $(srcdir)/BlockDavidsonEpetraEx.cpp

BlockDavidsonEpetraExGen_SOURCES = $(srcdir)/BlockDavidsonEpetraExGen.cpp

#np# Specify locations of include files within Trilinos, but outside np
#np# Add packages that np is dependent on; remove epetra if np does not need it

AM_CPPFLAGS = \
		-I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../belos/src -I$(top_builddir)/../belos/src \
		-I$(top_srcdir)/../ifpack/src -I$(top_builddir)/../ifpack/src \
		-I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src

#np# - List libraries that np is dependent on.  The DEPENDENCIES primary is
#np# optional.
BlockDavidsonEpetraEx_DEPENDENCIES = \
        $(top_builddir)/src/libanasazi.a \
        $(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a 

BlockDavidsonEpetraEx_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/util/ModeLaplace

#np# List libraries that np is dependent on.
BlockDavidsonEpetraEx_LDADD = \
        $(top_builddir)/src/libanasazi.a \
        $(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

BlockDavidsonEpetraExGen_DEPENDENCIES = \
        $(top_builddir)/src/libanasazi.a \
        $(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

BlockDavidsonEpetraExGen_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/util/ModeLaplace

#np# List libraries that np is dependent on.
BlockDavidsonEpetraExGen_LDADD = \
        $(top_builddir)/src/libanasazi.a \
        $(top_builddir)/util/ModeLaplace/libModeLaplace.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../ifpack/src/libifpack.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose examples are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# example1/ rather than this file.
# SUBDIRS = example1
SUBDIRS = 
