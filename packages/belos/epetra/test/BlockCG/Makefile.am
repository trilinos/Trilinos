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
BELOS_BCG_IFPACK_EX = bl_pcg_hb 
endif
BELOS_BCG_HB_EX = bl_cg_hb resolve_cg_hb Dummy
endif

noinst_PROGRAMS = $(BELOS_BCG_IFPACK_EX) $(BELOS_BCG_HB_EX)

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(top_srcdir)/epetra/example/BlockCG/bcsstk14.hb bcsstk14.hb
	cp  $(top_srcdir)/epetra/example/BlockCG/bcsstk15.hb bcsstk15.hb

Dummy_SOURCES = Dummy

bl_pcg_hb_SOURCES = $(srcdir)/test_bl_pcg_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_cg_hb_SOURCES = $(srcdir)/test_bl_cg_hb.cpp $(srcdir)/createEpetraProblem.cpp

resolve_cg_hb_SOURCES = $(srcdir)/test_resolve_cg_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_pcg_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_cg_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

resolve_cg_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_pcg_hb_CXXFLAGS = -I$(top_srcdir)/src

bl_cg_hb_CXXFLAGS = -I$(top_srcdir)/src 

resolve_cg_hb_CXXFLAGS = -I$(top_srcdir)/src 

bl_pcg_hb_LDADD = $(EXPORT_LIBS)

bl_cg_hb_LDADD = $(EXPORT_LIBS)

resolve_cg_hb_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 

EXTRA_DIST = \
	$(srcdir)/createEpetraProblem.hpp
