#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = ImportExport_test

ImportExport_test_SOURCES = $(srcdir)/cxx_main.cpp

#ImportExport_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

ImportExport_test_CXXFLAGS = -I$(top_srcdir)/src

if BUILD_BLAS
REFERENCE_BLAS_LIBS=-L$(top_builddir)/../../external/blas -lblas
endif

ImportExport_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(REFERENCE_BLAS_LIBS) $(LIBS) $(FLIBS)
