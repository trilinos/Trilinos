#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = bl_pcg_hb

bl_pcg_hb_SOURCES = $(srcdir)/test_bl_pcg_hb.cpp

bl_pcg_hb_DEPENDENCIES = $(top_builddir)/test/libBelosTest.a

bl_pcg_hb_CXXFLAGS = -I$(top_srcdir)/test -I$(top_srcdir)/src

bl_pcg_hb_LDADD = $(top_builddir)/test/libBelosTest.a $(LAPACK_LIBS) \
		$(BLAS_LIBS) $(LIBS) $(FLIBS)


