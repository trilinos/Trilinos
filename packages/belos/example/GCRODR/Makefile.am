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
BELOS_GCRODR_IFPACK_EX = GCRODRPrecEpetraExFile 
endif
BELOS_GCRODR_NOPREC_EX = GCRODREpetraExFile Dummy 
endif
endif

noinst_PROGRAMS = $(BELOS_GCRODR_NOPREC_EX) $(BELOS_GCRODR_IFPACK_EX)

Dummy$(EXEEXT): $(srcdir)/Dummy
	if [ $(srcdir) != '.' ]; then cp  $(srcdir)/sherman5.hb sherman5.hb; fi
	cp  $(srcdir)/../BlockGmres/orsirr1.hb orsirr1.hb
	cp  $(srcdir)/../BlockGmres/fidap036.hb fidap036.hb

Dummy_SOURCES = Dummy

GCRODREpetraExFile_SOURCES = $(srcdir)/GCRODREpetraExFile.cpp 

GCRODRPrecEpetraExFile_SOURCES = $(srcdir)/GCRODRPrecEpetraExFile.cpp 

GCRODREpetraExFile_DEPENDENCIES = \
			$(top_builddir)/../epetraext/src/libepetraext.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

GCRODRPrecEpetraExFile_DEPENDENCIES = \
			$(top_builddir)/../ifpack/src/libifpack.a \
			$(top_builddir)/../epetraext/src/libepetraext.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

GCRODREpetraExFile_CXXFLAGS = -I$(top_srcdir)/src 

GCRODRPrecEpetraExFile_CXXFLAGS = -I$(top_srcdir)/src 

GCRODREpetraExFile_LDADD = $(EXPORT_LIBS)

GCRODRPrecEpetraExFile_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 

EXTRA_DIST = \
        $(srcdir)/sherman5.hb
