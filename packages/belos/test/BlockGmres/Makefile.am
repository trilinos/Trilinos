#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = bl_pgmres_hb bl_gmres_hb

bl_pgmres_hb_SOURCES = $(srcdir)/test_bl_pgmres_hb.cpp

bl_gmres_hb_SOURCES = $(srcdir)/test_bl_gmres_hb.cpp

bl_pgmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a

bl_gmres_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a

bl_pgmres_hb_CXXFLAGS = -I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../ifpack/src \
		-I$(top_builddir)/../ifpack/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src

bl_gmres_hb_CXXFLAGS = 	-I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src

bl_pgmres_hb_LDADD = $(top_builddir)/../triutils/src/libtriutils.a \
		$(top_builddir)/../ifpack/src/libifpack.a \
		$(top_builddir)/../epetra/src/libepetra.a \
		$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

bl_gmres_hb_LDADD = $(top_builddir)/../triutils/src/libtriutils.a \
		$(top_builddir)/../epetra/src/libepetra.a \
		$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


#
