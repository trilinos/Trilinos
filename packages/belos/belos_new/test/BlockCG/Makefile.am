#@HEADER
#@HEADER

EXEEXT = .exe

#if BUILD_TRIUTILS
if BUILD_IFPACK
BELOS_BCG_IFPACK_EX = bl_pcg_hb
endif
if BUILD_TSFCORE_EPETRA
BELOS_BCG_TSFCORE_EX = bl_cg_tsfcore
if BUILD_IFPACK
BELOS_BCG_TSFCORE_IFPACK_EX = bl_pcg_tsfcore
endif
endif
BELOS_BCG_HB_EX = bl_cg_hb
#endif

noinst_PROGRAMS = $(BELOS_BCG_HB_EX) $(BELOS_BCG_IFPACK_EX) $(BELOS_BCG_TSFCORE_EX) $(BELOS_BCG_TSFCORE_IFPACK_EX) Dummy

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(srcdir)/bcsstk14.hb  bcsstk14.hb
	cp  $(srcdir)/bcsstk15.hb  bcsstk15.hb

Dummy_SOURCES = Dummy

bl_pcg_hb_SOURCES = $(srcdir)/test_bl_pcg_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_cg_hb_SOURCES = $(srcdir)/test_bl_cg_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_pcg_tsfcore_SOURCES = $(srcdir)/test_bl_pcg_tsfcore.cpp $(srcdir)/createEpetraProblem.cpp

bl_cg_tsfcore_SOURCES = $(srcdir)/test_bl_cg_tsfcore.cpp $(srcdir)/createEpetraProblem.cpp

bl_pcg_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../ifpack/src/libifpack.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_cg_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

bl_pcg_tsfcore_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
                        $(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
                        $(top_builddir)/../TSFCore/src/libTSFCore.a \
                        $(top_builddir)/../ifpack/src/libifpack.a \
                        $(top_builddir)/../epetra/src/libepetra.a \
                        $(top_builddir)/../teuchos/src/libteuchos.a

bl_cg_tsfcore_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
                        $(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
                        $(top_builddir)/../TSFCore/src/libTSFCore.a \
                        $(top_builddir)/../epetra/src/libepetra.a \
                        $(top_builddir)/../teuchos/src/libteuchos.a

# Include paths

common_include_paths = -I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../ifpack/src \
		-I$(top_builddir)/../ifpack/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src \
		-I$(top_srcdir)/../teuchos/src \
		-I$(top_builddir)/../teuchos/src

if BUILD_TSFCORE_EPETRA
common_include_paths += -I$(top_srcdir)/../TSFCore/src/interfaces/Core \
		-I$(top_srcdir)/../TSFCore/src/utilities/Core \
		-I$(top_srcdir)/../TSFCore/src/adapters/epetra \
		-I$(top_srcdir)/../TSFCore/src/adapters/mpi-base \
		-I$(top_srcdir)/../TSFCore/src/ \
		-I$(top_builddir)/../TSFCore/src \
		-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
		-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
		-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
		-I$(top_srcdir)/../TSFCoreUtils/src \
		-I$(top_builddir)/../TSFCoreUtils/src
endif

bl_pcg_hb_CXXFLAGS = $(common_include_paths)

bl_cg_hb_CXXFLAGS = $(common_include_paths)

bl_pcg_tsfcore_CXXFLAGS = $(common_include_paths)

bl_cg_tsfcore_CXXFLAGS = $(common_include_paths)

# Libraries

common_libraries = $(top_builddir)/../triutils/src/libtriutils.a

if BUILD_TSFCORE_EPETRA
common_libraries +=  $(top_builddir)/../TSFCore/src/libTSFCore.a \
                     $(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a
endif

if BUILD_IFPACK
common_libraries += $(top_builddir)/../ifpack/src/libifpack.a
endif

common_libraries += $(top_builddir)/../epetra/src/libepetra.a \
                $(top_builddir)/../teuchos/src/libteuchos.a \
                $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

bl_pcg_hb_LDADD = $(common_libraries)

bl_cg_hb_LDADD = $(common_libraries)

bl_cg_tsfcore_LDADD = $(common_libraries)

bl_pcg_tsfcore_LDADD = $(common_libraries)

#

