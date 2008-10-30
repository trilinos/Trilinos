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
BELOS_BCG_HB_EX = Tpetra_bl_cg_hb Dummy
endif

if BUILD_TRIUTILS
if BUILD_TEUCHOS_COMPLEX
BELOS_BCG_COMPLEX_TEST = Tpetra_bl_cg_complex_hb ComplexDummy
endif
endif

noinst_PROGRAMS = $(BELOS_BCG_HB_EX) $(BELOS_BCG_COMPLEX_TEST)

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(top_srcdir)/example/BlockCG/bcsstk14.hb bcsstk14.hb
	cp  $(top_srcdir)/example/BlockCG/bcsstk15.hb bcsstk15.hb

Dummy_SOURCES = Dummy

ComplexDummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(top_srcdir)/test/MVOPTester/mhd1280b.cua mhd1280b.cua

ComplexDummy_SOURCES = Dummy

Tpetra_bl_cg_hb_SOURCES = $(srcdir)/test_bl_cg_hb.cpp 

#resolve_cg_hb_SOURCES = $(srcdir)/test_resolve_cg_hb.cpp $(srcdir)/createTpetraProblem.cpp

Tpetra_bl_cg_complex_hb_SOURCES = $(srcdir)/test_bl_cg_complex_hb.cpp 

Tpetra_bl_cg_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../tpetra/src/libtpetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

#resolve_cg_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
#			$(top_builddir)/../tpetra/src/libtpetra.a \
#			$(top_builddir)/../teuchos/src/libteuchos.a

Tpetra_bl_cg_complex_hb_DEPENDENCIES = $(top_builddir)/src/libbelos.a \
			$(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../tpetra/src/libtpetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

Tpetra_bl_cg_hb_CXXFLAGS = -I$(top_srcdir)/src 

#resolve_cg_hb_CXXFLAGS = -I$(top_srcdir)/src 

Tpetra_bl_cg_complex_hb_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/test/MVOPTester

Tpetra_bl_cg_hb_LDADD = $(EXPORT_LIBS)

#resolve_cg_hb_LDADD = $(EXPORT_LIBS)

Tpetra_bl_cg_complex_hb_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 

EXTRA_DIST = # $(srcdir)/createTpetraProblem.hpp
