#@HEADER
#@HEADER

EXEEXT = .exe

if BUILD_TRIUTILS
if BUILD_IFPACK
BELOS_BGMRES_IFPACK_EX = bl_pgmres_hb bl_fgmres_hb
endif
BELOS_BGMRES_HB_EX = bl_gmres_hb
endif

noinst_PROGRAMS = $(BELOS_BGMRES_HB_EX) $(BELOS_BGMRES_IFPACK_EX) Dummy

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(srcdir)/orsirr1.hb orsirr1.hb
	cp  $(srcdir)/fidap036.hb fidap036.hb

Dummy_SOURCES = Dummy

bl_pgmres_hb_SOURCES = $(srcdir)/test_bl_pgmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_fgmres_hb_SOURCES = $(srcdir)/test_bl_fgmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

bl_gmres_hb_SOURCES = $(srcdir)/test_bl_gmres_hb.cpp $(srcdir)/createEpetraProblem.cpp

# Dependancies

bl_pgmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../ifpack/src/libifpack.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a			

bl_fgmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../ifpack/src/libifpack.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a			

bl_gmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
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

bl_pgmres_hb_CXXFLAGS = $(common_include_paths)

bl_fgmres_hb_CXXFLAGS = $(common_include_paths)

bl_gmres_hb_CXXFLAGS = $(common_include_paths)

# Libraries

common_libraries = $(top_builddir)/../triutils/src/libtriutils.a

if BUILD_IFPACK
common_libraries += $(top_builddir)/../ifpack/src/libifpack.a
endif

common_libraries += $(top_builddir)/../epetra/src/libepetra.a \
		$(top_builddir)/../teuchos/src/libteuchos.a \
		$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

bl_pgmres_hb_LDADD = $(common_libraries)

bl_fgmres_hb_LDADD = $(common_libraries)

bl_gmres_hb_LDADD = $(common_libraries)

#
