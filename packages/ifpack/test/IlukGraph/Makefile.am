#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = IlukGraph_test

IlukGraph_test_SOURCES = $(srcdir)/cxx_main.cpp

IlukGraph_test_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libifpack.a

IlukGraph_test_CXXFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src

IlukGraph_test_LDADD = -L$(top_builddir)/src -lifpack \
	-L$(top_builddir)/../epetra/src -lepetra $(FLIBS)


