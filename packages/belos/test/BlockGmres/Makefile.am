#@HEADER
#@HEADER

include $(top_builddir)/Makefile.export.belos

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(BELOS_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(BELOS_INCLUDES))
else
EXPORT_LIBS = $(BELOS_LIBS)
EXPORT_INC_PATH = $(BELOS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

#if BUILD_TRIUTILS
if BUILD_IFPACK
BELOS_BGMRES_IFPACK_EX = bl_pgmres_hb bl_fgmres_hb
endif
BELOS_BGMRES_EX = bl_gmres_hb bl_gmres_diag
#endif

noinst_PROGRAMS = $(BELOS_BGMRES_EX) $(BELOS_BGMRES_IFPACK_EX) Dummy

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(srcdir)/orsirr1.hb orsirr1.hb
	cp  $(srcdir)/fidap036.hb fidap036.hb

Dummy_SOURCES = Dummy

bl_pgmres_hb_SOURCES = $(srcdir)/test_bl_pgmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_fgmres_hb_SOURCES = $(srcdir)/test_bl_fgmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_gmres_hb_SOURCES = $(srcdir)/test_bl_gmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_gmres_diag_SOURCES = $(srcdir)/test_bl_gmres_diag.cpp

bl_pgmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a			

bl_fgmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a			

bl_gmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_gmres_diag_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_pgmres_hb_CXXFLAGS = -I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../ifpack/src \
		-I$(top_builddir)/../ifpack/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src \
		-I$(top_srcdir)/../teuchos/src \
		-I$(top_builddir)/../teuchos/src

bl_fgmres_hb_CXXFLAGS = -I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../ifpack/src \
		-I$(top_builddir)/../ifpack/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src \
		-I$(top_srcdir)/../teuchos/src \
		-I$(top_builddir)/../teuchos/src

bl_gmres_hb_CXXFLAGS = 	-I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src \
		-I$(top_srcdir)/../teuchos/src \
		-I$(top_builddir)/../teuchos/src

bl_gmres_diag_CXXFLAGS = -I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src \
		-I$(top_srcdir)/../teuchos/src \
		-I$(top_builddir)/../teuchos/src

bl_pgmres_hb_LDADD = $(EXPORT_LIBS)

bl_fgmres_hb_LDADD = $(EXPORT_LIBS)

bl_gmres_hb_LDADD = $(EXPORT_LIBS)

bl_gmres_diag_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 
#
