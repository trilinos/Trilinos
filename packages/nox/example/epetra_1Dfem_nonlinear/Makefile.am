## $Id$ 
## $Source$ 
#@HEADER
#@HEADER

if BUILD_NOXEPETRA
EXAMPLE = example
endif

noinst_PROGRAMS = $(EXAMPLE)

example_SOURCES = $(srcdir)/Example.C \
	$(srcdir)/FiniteElementProblem.C $(srcdir)/FiniteElementProblem.H \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H

example_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src/libnoxepetra.a

example_CXXFLAGS = $(TRILINOS_CXXFLAGS) \
	-I$(top_srcdir)/src

example_LDADD = -L$(top_builddir)/src -lnox -lnoxepetra 


