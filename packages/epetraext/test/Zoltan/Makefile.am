#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = Zoltan_test

Zoltan_test_SOURCES = \
  $(srcdir)/cxx_main.cpp \
  $(top_srcdir)/test/epetra_test_err.h

Zoltan_test_DEPENDENCIES = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/../epetra/src/libepetra.a

Zoltan_test_CXXFLAGS = -I$(top_srcdir)/src/transform \
                       -I$(top_builddir)/src \
                       -I$(top_srcdir)/../epetra/src \
                       -I$(top_builddir)/../epetra/src

Zoltan_test_LDADD = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/src/libepetra.a \
  $(BLAS_LIBS) \
  $(LIBS) \
  $(FLIBS)


