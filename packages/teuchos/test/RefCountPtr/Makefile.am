#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = RefCountPtr_test

RefCountPtr_test_SOURCES = $(srcdir)/cxx_main.cpp

RefCountPtr_test_DEPENDENCIES = $(top_builddir)/src/libteuchos.a

RefCountPtr_test_CXXFLAGS = -I$(top_srcdir)/src

RefCountPtr_test_LDADD = -L$(top_builddir)/src -lteuchos $(FLIBS)


