#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = SingletonFilter_test

SingletonFilter_test_SOURCES = \
  $(srcdir)/cxx_main.cpp \
  $(top_srcdir)/test/epetra_test_err.h

SingletonFilter_test_DEPENDENCIES = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/../epetra/src/libepetra.a

SingletonFilter_test_CXXFLAGS = -I$(top_srcdir)/src/transform \
                                -I$(top_builddir)/src \
                                -I$(top_srcdir)/../epetra/src \
                                -I$(top_builddir)/../epetra/src

SingletonFilter_test_LDADD = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/../epetra/src/libepetra.a \
  $(BLAS_LIBS) \
  $(LIBS) \
  $(FLIBS)


