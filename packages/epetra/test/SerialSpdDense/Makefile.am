#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = SerialSpdDense_test

SerialSpdDense_test_SOURCES = $(srcdir)/cxx_main.cpp

SerialSpdDense_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

SerialSpdDense_test_CXXFLAGS = -I$(top_srcdir)/src

SerialSpdDense_test_LDADD = -L$(top_builddir)/src -lepetra $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
