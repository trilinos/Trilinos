#@HEADER
#@HEADER

EXEEXT = .exe

#if BUILD_TRIUTILS
if BUILD_IFPACK
BELOS_BGMRES_IFPACK_EX = bl_pgmres_hb
endif
if BUILD_TSFCORE_EPETRA
BELOS_BGMRES_TSFCORE_EX = bl_gmres_tsfcore
endif
BELOS_BGMRES_HB_EX = bl_gmres_hb 
#endif

noinst_PROGRAMS = $(BELOS_BGMRES_HB_EX) $(BELOS_BGMRES_IFPACK_EX) $(BELOS_BGMRES_TSFCORE_EX) Dummy

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(srcdir)/orsirr1.hb orsirr1.hb
	cp  $(srcdir)/fidap036.hb fidap036.hb

Dummy_SOURCES = Dummy

bl_pgmres_hb_SOURCES = $(srcdir)/test_bl_pgmres_hb.cpp

bl_gmres_hb_SOURCES = $(srcdir)/test_bl_gmres_hb.cpp

bl_gmres_tsfcore_SOURCES = $(srcdir)/test_bl_gmres_tsfcore.cpp

bl_pgmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a			

bl_gmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_gmres_tsfcore_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
			$(top_builddir)/../TSFCore/src/libTSFCore.a \
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

bl_gmres_hb_CXXFLAGS = 	-I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src \
		-I$(top_srcdir)/../teuchos/src \
		-I$(top_builddir)/../teuchos/src

bl_gmres_tsfcore_CXXFLAGS = 	-I$(top_srcdir)/src \
		-I$(top_srcdir)/../TSFCore/src/interfaces/Core \
		-I$(top_srcdir)/../TSFCore/src/utilities/Core \
		-I$(top_srcdir)/../TSFCore/src/adapters/epetra \
		-I$(top_srcdir)/../TSFCore/src/adapters/mpi-base \
		-I$(top_srcdir)/../TSFCore/src/ \
		-I$(top_builddir)/../TSFCore/src \
		-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
		-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
		-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
		-I$(top_srcdir)/../TSFCoreUtils/src \
		-I$(top_builddir)/../TSFCoreUtils/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src \
		-I$(top_srcdir)/../teuchos/src \
		-I$(top_builddir)/../teuchos/src

bl_pgmres_hb_LDADD = $(top_builddir)/../triutils/src/libtriutils.a \
		$(top_builddir)/../ifpack/src/libifpack.a \
		$(top_builddir)/../epetra/src/libepetra.a \
		$(top_builddir)/../teuchos/src/libteuchos.a \
		$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

bl_gmres_hb_LDADD = $(top_builddir)/../triutils/src/libtriutils.a \
		$(top_builddir)/../epetra/src/libepetra.a \
		$(top_builddir)/../teuchos/src/libteuchos.a \
		$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

bl_gmres_tsfcore_LDADD = $(top_builddir)/../triutils/src/libtriutils.a \
		$(top_builddir)/../TSFCore/src/libTSFCore.a \
		$(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
		$(top_builddir)/../epetra/src/libepetra.a \
		$(top_builddir)/../teuchos/src/libteuchos.a \
		$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

#
