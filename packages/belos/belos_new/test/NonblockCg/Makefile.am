#@HEADER
#@HEADER

EXEEXT = .exe

BELOS_NBCG_TSFCORE_SIMPLE_EX = nbl_pcg_tsfcore_simple

#if BUILD_TRIUTILS
#if BUILD_TSFCORE_EPETRA
#if BUILD_IFPACK
#BELOS_NBCG_TSFCORE_IFPACK_EX = nbl_pcg_tsfcore
#endif
#endif
#endif


Dummy$(EXEEXT):
	cp  $(srcdir)/run_simple_tests.pl ./run_simple_tests.pl

#noinst_PROGRAMS = $(BELOS_NBCG_TSFCORE_IFPACK_EX) $(BELOS_NBCG_TSFCORE_SIMPLE_EX)
noinst_PROGRAMS = $(BELOS_NBCG_TSFCORE_SIMPLE_EX) Dummy

#nbl_pcg_tsfcore_SOURCES = $(srcdir)/test_nbl_pcg_tsfcore.cpp $(srcdir)/../BlockCG/createEpetraProblem.cpp

nbl_pcg_tsfcore_simple_SOURCES = $(srcdir)/test_nbl_pcg_tsfcore_simple.cpp

#nbl_pcg_tsfcore_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
#                        $(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
#                        $(top_builddir)/../TSFCore/src/libTSFCore.a \
#                        $(top_builddir)/../ifpack/src/libifpack.a \
#                        $(top_builddir)/../epetra/src/libepetra.a \
#                        $(top_builddir)/../teuchos/src/libteuchos.a


nbl_pcg_tsfcore_simple_DEPENDENCIES = $(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a \
                        $(top_builddir)/../TSFCore/src/libTSFCore.a \
                        $(top_builddir)/../teuchos/src/libteuchos.a

# Include paths

common_include_paths = -I$(top_srcdir)/src \
	-I$(top_srcdir)/src/interfaces \
	-I$(top_srcdir)/src/implementations \
	-I$(srcdir)/../BlockCG \
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
		-I$(top_srcdir)/../TSFCore/src/interfaces/Solvers \
		-I$(top_srcdir)/../TSFCore/src/utilities/Core \
		-I$(top_srcdir)/../TSFCore/src/utilities/Solvers \
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

#nbl_pcg_tsfcore_CXXFLAGS = $(common_include_paths)

nbl_pcg_tsfcore_simple_CXXFLAGS = $(common_include_paths)

# Libraries

common_libraries =  $(top_builddir)/../TSFCore/src/libTSFCore.a \
                     $(top_builddir)/../TSFCoreUtils/src/libTSFCoreUtils.a

#if BUILD_IFPACK
#common_libraries += $(top_builddir)/../ifpack/src/libifpack.a
#endif

common_libraries += $(top_builddir)/../epetra/src/libepetra.a \
                $(top_builddir)/../teuchos/src/libteuchos.a \
                $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

#nbl_pcg_tsfcore_LDADD = $(common_libraries)

nbl_pcg_tsfcore_simple_LDADD = $(common_libraries)

#

