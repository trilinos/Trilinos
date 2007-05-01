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
BELOS_BCG_IFPACK_EX = BlockPrecCGEpetraExFile
endif
BELOS_BCG_NOPREC_EX = BlockCGEpetraExFile Dummy 
endif
endif

#noinst_PROGRAMS = $(BELOS_BGMRES_EX) $(BELOS_BGMRES_IFPACK_EX) $(BELOS_BGMRES_NOPREC_EX) $(BELOS_BGMRES_COMPLEX_EX)
noinst_PROGRAMS = $(BELOS_BCG_NOPREC_EX) $(BELOS_BCG_IFPACK_EX)

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(srcdir)/bcsstk14.hb bcsstk14.hb
	cp  $(srcdir)/bcsstk15.hb bcsstk15.hb

Dummy_SOURCES = Dummy

BlockPrecCGEpetraExFile_SOURCES = $(srcdir)/BlockPrecCGEpetraExFile.cpp 

BlockCGEpetraExFile_SOURCES = $(srcdir)/BlockCGEpetraExFile.cpp 

BlockCGEpetraExFile_DEPENDENCIES = \
			$(top_builddir)/../epetraext/src/libepetraext.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

BlockPrecCGEpetraExFile_DEPENDENCIES = \
			$(top_builddir)/../ifpack/src/libifpack.a \
			$(top_builddir)/../epetraext/src/libepetraext.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

BlockCGEpetraExFile_CXXFLAGS = -I$(top_srcdir)/src 

BlockPrecCGEpetraExFile_CXXFLAGS = -I$(top_srcdir)/src 

BlockCGEpetraExFile_LDADD = $(EXPORT_LIBS)

BlockPrecCGEpetraExFile_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 
#
