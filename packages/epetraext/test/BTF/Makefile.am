#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = BTF_test

BTF_test_SOURCES = \
  $(srcdir)/cxx_main.cpp \
  $(top_srcdir)/test/epetra_test_err.h

BTF_test_DEPENDENCIES = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/../epetra/src/libepetra.a

BTF_test_CXXFLAGS = -I$(top_builddir)/src \
                    -I$(top_srcdir)/src/transform \
                    -I$(top_builddir)/../epetra/src \
                    -I$(top_srcdir)/../epetra/src

BTF_test_LDADD = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/../epetra/src/libepetra.a \
  $(BLAS_LIBS) \
  $(LIBS) \
  $(FLIBS)


