#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = VbrMatrix_test

VbrMatrix_test_SOURCES = $(srcdir)/cxx_main.cpp

VbrMatrix_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

VbrMatrix_test_CXXFLAGS = -I$(top_srcdir)/src

VbrMatrix_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
