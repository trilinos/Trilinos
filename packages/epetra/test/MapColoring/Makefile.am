#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = MapColoring_test

MapColoring_test_SOURCES = $(srcdir)/cxx_main.cpp

MapColoring_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

MapColoring_test_CXXFLAGS = -I$(top_srcdir)/src

MapColoring_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)


