#@HEADER
#@HEADER

EXAMPLE = power_method_example

noinst_PROGRAMS = $(EXAMPLE)

power_method_example_SOURCES = $(srcdir)/cxx_main.cpp

#power_method_example_DEPENDENCIES = $(top_builddir)/src/libepetra.a

power_method_example_CXXFLAGS = -I$(top_srcdir)/src

if BUILD_BLAS
REFERENCE_BLAS_LIBS=-L$(top_builddir)/../../external/blas -lblas
endif


power_method_example_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(REFERENCE_BLAS_LIBS) $(LIBS) $(FLIBS)
