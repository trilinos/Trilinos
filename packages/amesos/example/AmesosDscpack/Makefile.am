#@HEADER
#@HEADER

EXAMPLE = example

noinst_PROGRAMS = $(EXAMPLE)

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src

example_SOURCES = $(srcdir)/CallAmesosDscpack.cpp \
	$(srcdir)/CreateTridi.cpp $(srcdir)/CreateTridi.h

example_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/amd/src/libamd.a \
	$(top_builddir)/../triutils/src/libtriutils.a 

example_CXXFLAGS = -I$(top_srcdir)/src

example_LDADD = $(top_builddir)/src/libamesos.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(AMESOS_LIBS) \
        $(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS) 

