#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = bl_pcg_hb bl_cg_hb

bl_pcg_hb_SOURCES = $(srcdir)/test_bl_pcg_hb.cpp 

bl_cg_hb_SOURCES = $(srcdir)/test_bl_cg_hb.cpp 

bl_pcg_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a

bl_cg_hb_DEPENDENCIES = $(top_builddir)/../triutils/src/libtriutils.a \
			$(top_builddir)/../epetra/src/libepetra.a

bl_pcg_hb_CXXFLAGS = -I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../ifpack/src \
		-I$(top_builddir)/../ifpack/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src

bl_cg_hb_CXXFLAGS = -I$(top_srcdir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src

bl_pcg_hb_LDADD = $(top_builddir)/../epetra/src/libepetra.a \
		$(top_builddir)/../triutils/src/libtriutils.a \
		$(top_builddir)/../ifpack/src/libifpack.a $(LAPACK_LIBS) \
		$(BLAS_LIBS) $(LIBS) $(FLIBS)

bl_cg_hb_LDADD = $(top_builddir)/../epetra/src/libepetra.a \
		$(top_builddir)/../triutils/src/libtriutils.a \
		$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


