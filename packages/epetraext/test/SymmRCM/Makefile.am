#@HEADER
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = SymmRCM_test

SymmRCM_test_SOURCES = \
  $(srcdir)/cxx_main.cpp \
  $(top_srcdir)/test/epetra_test_err.h

SymmRCM_test_DEPENDENCIES = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/../epetra/src/libepetra.a

SymmRCM_test_CXXFLAGS = -I$(top_srcdir)/src/transform \
                        -I$(top_builddir)/src \
                        -I$(top_srcdir)/../epetra/src \
                        -I$(top_builddir)/../epetra/src

SymmRCM_test_LDADD = \
  $(top_builddir)/src/transform/libepetraext_trans.a \
  $(top_builddir)/../epetra/src/libepetra.a \
  $(BLAS_LIBS) \
  $(LIBS) \
  $(FLIBS)


