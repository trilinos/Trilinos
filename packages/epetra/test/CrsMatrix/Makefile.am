#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = CrsMatrix_test

CrsMatrix_test_SOURCES = $(srcdir)/cxx_main.cpp

CrsMatrix_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

CrsMatrix_test_CXXFLAGS = -I$(top_srcdir)/src

CrsMatrix_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
