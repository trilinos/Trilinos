#@HEADER
#@HEADER

EXAMPLE = example

noinst_PROGRAMS = $(EXAMPLE)

example_SOURCES = $(srcdir)/cxx_main.cpp

#example_DEPENDENCIES = $(top_builddir)/src/libepetra.a

example_CXXFLAGS = -I$(top_srcdir)/src

if BUILD_BLAS
REFERENCE_BLAS_LIBS=-L$(top_builddir)/../../external/blas -lblas
endif

example_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(REFERENCE_BLAS_LIBS) $(LIBS) $(FLIBS)
