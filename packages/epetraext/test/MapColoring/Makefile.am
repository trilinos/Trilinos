#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = MapColoring_test

MapColoring_test_SOURCES = \
  $(srcdir)/cxx_main.cpp \
  $(top_srcdir)/test/epetra_test_err.h

MapColoring_test_DEPENDENCIES = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/../epetra/src/libepetra.a

MapColoring_test_CXXFLAGS = -I$(top_srcdir)/src/transform \
                            -I$(top_builddir)/src \
                            -I$(top_srcdir)/../epetra/src \
                            -I$(top_builddir)/../epetra/src

MapColoring_test_LDADD = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/../epetra/src/libepetra.a \
  $(BLAS_LIBS) \
  $(LIBS) \
  $(FLIBS)

