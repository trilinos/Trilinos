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
BELOS_BGMRES_COMPLEX_EX = bl_gmres_complex_hb ComplexDummy
endif
endif


noinst_PROGRAMS = $(BELOS_BGMRES_COMPLEX_EX)

ComplexDummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(top_srcdir)/test/MVOPTester/mhd1280b.cua mhd1280b.cua
	cp  $(top_srcdir)/test/MVOPTester/mhd1280a.cua mhd1280a.cua

ComplexDummy_SOURCES = Dummy

bl_gmres_complex_hb_SOURCES = $(srcdir)/test_bl_gmres_complex_hb.cpp

bl_gmres_complex_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_gmres_complex_hb_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/test/MVOPTester

bl_gmres_complex_hb_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 
