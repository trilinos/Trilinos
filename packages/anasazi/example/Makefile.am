#@HEADER
#@HEADER

EXEEXT = .exe

if BUILD_BELOS
if BUILD_IFPACK
ANASAZI_BIRAM_EX_SI_IFPACK = BlockArnoldiPetraExGen
endif
endif

#if BUILD_TRIUTILS
ANASAZI_BIRAM_EX_HB = BlockArnoldiPetraExHb
#endif 

noinst_PROGRAMS = $(ANASAZI_BIRAM_EX_SI_IFPACK) $(ANASAZI_BIRAM_EX_HB) BlockArnoldiPetraEx BlockArnoldiPetraExSVD

BlockArnoldiPetraEx_SOURCES = $(srcdir)/BlockArnoldiPetraEx.cpp

BlockArnoldiPetraExSVD_SOURCES = $(srcdir)/BlockArnoldiPetraExSVD.cpp

BlockArnoldiPetraExHb_SOURCES = $(srcdir)/BlockArnoldiPetraExHb.cpp

BlockArnoldiPetraExGen_SOURCES = $(srcdir)/BlockArnoldiPetraExGen.cpp

#np# Specify locations of include files within Trilinos, but outside np
#np# Add packages that np is dependent on; remove epetra if np does not need it

AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../belos/src -I$(top_builddir)/../belos/src \
		-I$(top_srcdir)/../ifpack/src -I$(top_builddir)/../ifpack/src \
		-I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src

#np# - List libraries that np is dependent on.  The DEPENDENCIES primary is
#np# optional.
BlockArnoldiPetraEx_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a 

BlockArnoldiPetraEx_CXXFLAGS = -I$(top_srcdir)/src

BlockArnoldiPetraEx_LDADD = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

BlockArnoldiPetraExSVD_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a 

BlockArnoldiPetraExSVD_CXXFLAGS = -I$(top_srcdir)/src 

BlockArnoldiPetraExSVD_LDADD = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

BlockArnoldiPetraExHb_DEPENDENCIES = \
	$(top_builddir)/../triutils/src/libtriutils.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

BlockArnoldiPetraExHb_CXXFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../triutils/src

#np# List libraries that np is dependent on.
BlockArnoldiPetraExHb_LDADD = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../triutils/src/libtriutils.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

BlockArnoldiPetraExGen_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

BlockArnoldiPetraExGen_CXXFLAGS = -I$(top_srcdir)/src

#np# List libraries that np is dependent on.
BlockArnoldiPetraExGen_LDADD = \
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
