#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = CrsRiluk_test

CrsRiluk_test_SOURCES = $(srcdir)/cxx_main.cpp

#CrsRiluk_test_DEPENDENCIES = \
#	$(top_builddir)/../epetra/src/libepetra.a \
#	$(top_builddir)/src/libifpack.a

CrsRiluk_test_CXXFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src

if BUILD_BLAS
REFERENCE_BLAS_LIBS=-L$(top_builddir)/../../external/blas -lblas
endif

CrsRiluk_test_LDADD = -L$(top_builddir)/src -lifpack \
	-L$(top_builddir)/../epetra/src -lepetra \
	$(BLAS_LIBS) $(REFERENCE_BLAS_LIBS) $(FLIBS)


