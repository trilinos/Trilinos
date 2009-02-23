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
BELOS_BCG_HB_EX = Tpetra_bl_gmres_hb_df Dummy
endif

noinst_PROGRAMS = $(BELOS_BCG_HB_EX) $(BELOS_BCG_COMPLEX_TEST)

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(top_srcdir)/tpetra/test/BlockCG/bcsstk17.rsa bcsstk17.rsa
	cp  $(top_srcdir)/tpetra/test/BlockCG/bcsstk18.rsa bcsstk18.rsa

Dummy_SOURCES = Dummy

ComplexDummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(top_srcdir)/test/MVOPTester/mhd1280b.cua mhd1280b.cua

Tpetra_bl_gmres_hb_df_SOURCES = $(srcdir)/test_bl_gmres_hb_df.cpp 

Tpetra_bl_gmres_hb_df_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../tpetra/src/libtpetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

Tpetra_bl_gmres_hb_df_CXXFLAGS = -I$(top_srcdir)/src 

Tpetra_bl_gmres_hb_df_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 

EXTRA_DIST = # $(srcdir)/createTpetraProblem.hpp
