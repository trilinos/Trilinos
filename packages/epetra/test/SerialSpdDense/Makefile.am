#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = SerialSpdDense_test

SerialSpdDense_test_SOURCES = $(srcdir)/cxx_main.cpp

#SerialSpdDense_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

SerialSpdDense_test_CXXFLAGS = -I$(top_srcdir)/src

if BUILD_BLAS
REFERENCE_BLAS_LIBS=-L$(top_builddir)/../../external/blas -lblas
endif

if BUILD_LAPACK
REFERENCE_LAPACK_LIBS=-L$(top_builddir)/../../external/lapack -llapack
endif

SerialSpdDense_test_LDADD = -L$(top_builddir)/src -lepetra $(LAPACK_LIBS) $(REFERENCE_LAPACK_LIBS) $(BLAS_LIBS) $(REFERENCE_BLAS_LIBS) $(LIBS) $(FLIBS)
