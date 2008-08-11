#@HEADER
#@HEADER

include $(top_builddir)/Makefile.export.belos

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(BELOS_TEST_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(BELOS_TEST_INCLUDES))
else
EXPORT_LIBS = $(BELOS_TEST_LIBS)
EXPORT_INC_PATH = $(BELOS_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

if BUILD_TRIUTILS
if BUILD_EPETRAEXT
if BUILD_IFPACK
BELOS_BGMRES_IFPACK_EX = BlockPrecGmresEpetraExFile BlockFlexGmresEpetraExFile PseudoBlockPrecGmresEpetraExFile
endif
BELOS_BGMRES_NOPREC_EX = BlockGmresEpetraExFile BlockGmresPolyEpetraExFile PseudoBlockGmresEpetraExFile Dummy 
endif
endif

noinst_PROGRAMS = $(BELOS_BGMRES_NOPREC_EX) $(BELOS_BGMRES_IFPACK_EX)

Dummy$(EXEEXT): $(srcdir)/Dummy
	if [ $(srcdir) != '.' ]; then cp  $(srcdir)/orsirr1.hb orsirr1.hb; fi
	if [ $(srcdir) != '.' ]; then cp  $(srcdir)/fidap036.hb fidap036.hb; fi

Dummy_SOURCES = Dummy

BlockPrecGmresEpetraExFile_SOURCES = $(srcdir)/BlockPrecGmresEpetraExFile.cpp 

BlockFlexGmresEpetraExFile_SOURCES = $(srcdir)/BlockFlexGmresEpetraExFile.cpp 

BlockGmresEpetraExFile_SOURCES = $(srcdir)/BlockGmresEpetraExFile.cpp 

BlockGmresPolyEpetraExFile_SOURCES = $(srcdir)/BlockGmresPolyEpetraExFile.cpp 

PseudoBlockPrecGmresEpetraExFile_SOURCES = $(srcdir)/PseudoBlockPrecGmresEpetraExFile.cpp 

PseudoBlockGmresEpetraExFile_SOURCES = $(srcdir)/PseudoBlockGmresEpetraExFile.cpp 

BlockGmresEpetraExFile_DEPENDENCIES = \
			$(top_builddir)/../epetraext/src/libepetraext.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

BlockGmresPolyEpetraExFile_DEPENDENCIES = \
			$(top_builddir)/../epetraext/src/libepetraext.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

BlockPrecGmresEpetraExFile_DEPENDENCIES = \
			$(top_builddir)/../ifpack/src/libifpack.a \
			$(top_builddir)/../epetraext/src/libepetraext.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

BlockFlexGmresEpetraExFile_DEPENDENCIES = \
			$(top_builddir)/../ifpack/src/libifpack.a \
			$(top_builddir)/../epetraext/src/libepetraext.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

PseudoBlockGmresEpetraExFile_DEPENDENCIES = \
			$(top_builddir)/../epetraext/src/libepetraext.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

PseudoBlockPrecGmresEpetraExFile_DEPENDENCIES = \
			$(top_builddir)/../ifpack/src/libifpack.a \
			$(top_builddir)/../epetraext/src/libepetraext.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

BlockGmresEpetraExFile_CXXFLAGS = -I$(top_srcdir)/src 

BlockGmresPolyEpetraExFile_CXXFLAGS = -I$(top_srcdir)/src 

BlockPrecGmresEpetraExFile_CXXFLAGS = -I$(top_srcdir)/src 

BlockFlexGmresEpetraExFile_CXXFLAGS = -I$(top_srcdir)/src 

PseudoBlockGmresEpetraExFile_CXXFLAGS = -I$(top_srcdir)/src 

PseudoBlockPrecGmresEpetraExFile_CXXFLAGS = -I$(top_srcdir)/src 

BlockGmresEpetraExFile_LDADD = $(EXPORT_LIBS)

BlockGmresPolyEpetraExFile_LDADD = $(EXPORT_LIBS)

BlockPrecGmresEpetraExFile_LDADD = $(EXPORT_LIBS)

BlockFlexGmresEpetraExFile_LDADD = $(EXPORT_LIBS)

PseudoBlockGmresEpetraExFile_LDADD = $(EXPORT_LIBS)

PseudoBlockPrecGmresEpetraExFile_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 
#
EXTRA_DIST = \
        $(srcdir)/orsirr1.hb \
        $(srcdir)/fidap036.hb
