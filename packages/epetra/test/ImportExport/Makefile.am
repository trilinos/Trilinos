#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = ImportExport_test

ImportExport_test_SOURCES = $(srcdir)/cxx_main.cpp

ImportExport_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

ImportExport_test_CXXFLAGS = -I$(top_srcdir)/src

ImportExport_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
