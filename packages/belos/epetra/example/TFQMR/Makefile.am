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
if BUILD_IFPACK
BELOS_TFQMR_PREC_EX = TFQMREpetraExFile Dummy
endif
endif

noinst_PROGRAMS = $(BELOS_TFQMR_PREC_EX)

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(top_srcdir)/epetra/example/BlockGmres/orsirr1.hb orsirr1.hb 

Dummy_SOURCES = Dummy

TFQMREpetraExFile_SOURCES = $(srcdir)/TFQMREpetraExFile.cpp $(srcdir)/createEpetraProblem.cpp

TFQMREpetraExFile_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../ifpack/src/libifpack.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

TFQMREpetraExFile_CXXFLAGS = 	-I$(top_srcdir)/src 

TFQMREpetraExFile_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 
#
