#@HEADER
#@HEADER

EXAMPLE = power_method_example

noinst_PROGRAMS = $(EXAMPLE)

power_method_example_SOURCES = $(srcdir)/cxx_main.cpp

power_method_example_DEPENDENCIES = $(top_builddir)/src/libepetra.a

power_method_example_CXXFLAGS = -I$(top_srcdir)/src

power_method_example_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
