#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = IntSerialDense_test

IntSerialDense_test_SOURCES = $(srcdir)/cxx_main.cpp

IntSerialDense_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

IntSerialDense_test_CXXFLAGS = -I$(top_srcdir)/src

IntSerialDense_test_LDADD = -L$(top_builddir)/src -lepetra $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


