#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = BlockMap_test

BlockMap_test_SOURCES = $(srcdir)/cxx_main.cpp $(srcdir)/checkmap.cpp $(srcdir)/checkmap.h $(top_srcdir)/test/epetra_test_err.h

BlockMap_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

BlockMap_test_CXXFLAGS = -I$(top_srcdir)/src

BlockMap_test_LDADD = -L$(top_builddir)/src -lepetra $(FLIBS)


