#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = Map_test

Map_test_SOURCES = $(srcdir)/cxx_main.cpp $(srcdir)/checkmap.cpp $(srcdir)/checkmap.h

Map_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

Map_test_CXXFLAGS = -I$(top_srcdir)/src
#Unsure why this test has to link with BLAS for IBM
Map_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)


