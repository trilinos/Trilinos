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
BELOS_TFQMR_NOPREC_EX = tfqmr_hb Dummy 
if BUILD_IFPACK
BELOS_TFQMR_PREC_EX = ptfqmr_hb
endif
endif

noinst_PROGRAMS = $(BELOS_TFQMR_NOPREC_EX) $(BELOS_TFQMR_PREC_EX)

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(top_srcdir)/example/BlockGmres/orsirr1.hb orsirr1.hb 
	cp  $(top_srcdir)/example/BlockGmres/fidap036.hb fidap036.hb

Dummy_SOURCES = Dummy

tfqmr_hb_SOURCES = $(srcdir)/test_tfqmr_hb.cpp $(srcdir)/createEpetraProblem.cpp

tfqmr_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

tfqmr_hb_CXXFLAGS = 	-I$(top_srcdir)/src 

tfqmr_hb_LDADD = $(EXPORT_LIBS)

ptfqmr_hb_SOURCES = $(srcdir)/test_ptfqmr_hb.cpp $(srcdir)/createEpetraProblem.cpp

ptfqmr_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../ifpack/src/libifpack.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

ptfqmr_hb_CXXFLAGS = 	-I$(top_srcdir)/src 

ptfqmr_hb_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 
#
