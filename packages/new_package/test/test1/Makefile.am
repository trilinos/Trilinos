#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = Object_test

Object_test_SOURCES = $(srcdir)/cxx_main.cpp

Object_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

Object_test_CXXFLAGS = -I$(top_srcdir)/src

Object_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)


