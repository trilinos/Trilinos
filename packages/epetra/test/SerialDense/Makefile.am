#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = SerialDense_test

SerialDense_test_SOURCES = $(srcdir)/cxx_main.cpp

SerialDense_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

SerialDense_test_CXXFLAGS = -I$(top_srcdir)/src

SerialDense_test_LDADD = -L$(top_builddir)/src -lepetra $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


