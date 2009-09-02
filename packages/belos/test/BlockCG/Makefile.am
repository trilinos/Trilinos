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
if BUILD_TEUCHOS_COMPLEX
BELOS_BCG_COMPLEX_TEST = bl_cg_complex_hb ComplexDummy
endif
endif

noinst_PROGRAMS = $(BELOS_BCG_COMPLEX_TEST)

ComplexDummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(top_srcdir)/test/MVOPTester/mhd1280b.cua mhd1280b.cua

ComplexDummy_SOURCES = Dummy

bl_cg_complex_hb_SOURCES = $(srcdir)/test_bl_cg_complex_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_cg_complex_hb_DEPENDENCIES = $(top_builddir)/src/libbelos.a \
			$(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_cg_complex_hb_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/test/MVOPTester

bl_cg_complex_hb_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 
