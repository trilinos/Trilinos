#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = View_test

View_test_SOURCES = \
  $(srcdir)/cxx_main.cpp \
  $(top_srcdir)/test/epetra_test_err.h

View_test_DEPENDENCIES = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/../epetra/src/libepetra.a

View_test_CXXFLAGS = -I$(top_srcdir)/src/transform \
                     -I$(top_builddir)/src \
                     -I$(top_srcdir)/../epetra/src \
                     -I$(top_builddir)/../epetra/src

View_test_LDADD = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/../epetra/src/libepetra.a \
  $(BLAS_LIBS) \
  $(LIBS) \
  $(FLIBS)


