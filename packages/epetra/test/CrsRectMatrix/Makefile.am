#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = CrsRectMatrix_test

CrsRectMatrix_test_SOURCES = $(srcdir)/cxx_main.cpp

CrsRectMatrix_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

CrsRectMatrix_test_CXXFLAGS = -I$(top_srcdir)/src

CrsRectMatrix_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
