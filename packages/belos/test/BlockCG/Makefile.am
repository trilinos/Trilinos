#@HEADER
#@HEADER

EXEEXT = .exe

#if BUILD_TRIUTILS
if BUILD_IFPACK
BELOS_BCG_IFPACK_EX = bl_pcg_hb
endif
BELOS_BCG_HB_EX = bl_cg_hb
#endif

noinst_PROGRAMS = $(BELOS_BCG_HB_EX) $(BELOS_BCG_IFPACK_EX) Dummy

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(srcdir)/bcsstk14.hb  bcsstk14.hb
	cp  $(srcdir)/bcsstk15.hb  bcsstk15.hb

Dummy_SOURCES = Dummy

bl_pcg_hb_SOURCES = $(srcdir)/test_bl_pcg_hb.cpp 

bl_cg_hb_SOURCES = $(srcdir)/test_bl_cg_hb.cpp 

bl_pcg_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_cg_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_pcg_hb_CXXFLAGS = -I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../ifpack/src \
		-I$(top_builddir)/../ifpack/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src \
		-I$(top_srcdir)/../teuchos/src \
		-I$(top_builddir)/../teuchos/src

bl_cg_hb_CXXFLAGS = -I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src \
		-I$(top_srcdir)/../teuchos/src \
		-I$(top_builddir)/../teuchos/src

bl_pcg_hb_LDADD = $(top_builddir)/../epetra/src/libepetra.a \
		$(top_builddir)/../triutils/src/libtriutils.a \
		$(top_builddir)/../ifpack/src/libifpack.a  \
		$(top_builddir)/../teuchos/src/libteuchos.a \
		$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

bl_cg_hb_LDADD = $(top_builddir)/../epetra/src/libepetra.a \
		$(top_builddir)/../triutils/src/libtriutils.a \
		$(top_builddir)/../teuchos/src/libteuchos.a \
		$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


