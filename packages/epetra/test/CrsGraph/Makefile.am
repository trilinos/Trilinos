#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = CrsGraph_test

CrsGraph_test_SOURCES = $(srcdir)/cxx_main.cpp

CrsGraph_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

CrsGraph_test_CXXFLAGS = -I$(top_srcdir)/src

CrsGraph_test_LDADD = -L$(top_builddir)/src -lepetra $(FLIBS)


