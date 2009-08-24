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
BELOS_BGMRES_IFPACK_EX = bl_pgmres_hb pseudo_pgmres_hb bl_fgmres_hb
endif
BELOS_BGMRES_NOPREC_EX = bl_gmres_hb resolve_gmres_hb pseudo_gmres_hb Dummy 
endif
BELOS_BGMRES_EX = bl_gmres_diag

noinst_PROGRAMS = $(BELOS_BGMRES_EX) $(BELOS_BGMRES_IFPACK_EX) $(BELOS_BGMRES_NOPREC_EX) 

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(top_srcdir)/epetra/example/BlockGmres/orsirr1.hb orsirr1.hb
	cp  $(top_srcdir)/epetra/example/BlockGmres/fidap036.hb fidap036.hb

Dummy_SOURCES = Dummy


bl_pgmres_hb_SOURCES = $(srcdir)/test_bl_pgmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

pseudo_pgmres_hb_SOURCES = $(srcdir)/test_pseudo_pgmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_fgmres_hb_SOURCES = $(srcdir)/test_bl_fgmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_gmres_hb_SOURCES = $(srcdir)/test_bl_gmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

resolve_gmres_hb_SOURCES = $(srcdir)/test_resolve_gmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

pseudo_gmres_hb_SOURCES = $(srcdir)/test_pseudo_gmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_gmres_diag_SOURCES = $(srcdir)/test_bl_gmres_diag.cpp

bl_pgmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a			

pseudo_pgmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a			

bl_fgmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a			

bl_gmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

resolve_gmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

pseudo_gmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_gmres_diag_DEPENDENCIES = \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_pgmres_hb_CXXFLAGS = -I$(top_srcdir)/src

pseudo_pgmres_hb_CXXFLAGS = -I$(top_srcdir)/src

bl_fgmres_hb_CXXFLAGS = -I$(top_srcdir)/src 

bl_gmres_hb_CXXFLAGS = 	-I$(top_srcdir)/src 

resolve_gmres_hb_CXXFLAGS =  -I$(top_srcdir)/src 

pseudo_gmres_hb_CXXFLAGS = -I$(top_srcdir)/src 

bl_gmres_diag_CXXFLAGS = -I$(top_srcdir)/src

bl_pgmres_hb_LDADD = $(EXPORT_LIBS)

pseudo_pgmres_hb_LDADD = $(EXPORT_LIBS)

bl_fgmres_hb_LDADD = $(EXPORT_LIBS)

bl_gmres_hb_LDADD = $(EXPORT_LIBS)

resolve_gmres_hb_LDADD = $(EXPORT_LIBS)

pseudo_gmres_hb_LDADD = $(EXPORT_LIBS)

bl_gmres_diag_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 

EXTRA_DIST = \
	$(srcdir)/createEpetraProblem.hpp
#
