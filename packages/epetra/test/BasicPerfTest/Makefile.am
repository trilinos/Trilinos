#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = BasicPerfTest_test

BasicPerfTest_test_SOURCES = $(srcdir)/cxx_main.cpp

BasicPerfTest_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

BasicPerfTest_test_CXXFLAGS = -I$(top_srcdir)/src

BasicPerfTest_test_LDADD = $(top_builddir)/src/libepetra.a $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
