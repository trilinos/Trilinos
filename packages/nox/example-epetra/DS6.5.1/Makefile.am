## $Id$ 
## $Source$ 
#@HEADER
#@HEADER

if BUILD_NOXEPETRA
EXAMPLE = example
endif

noinst_PROGRAMS = $(EXAMPLE)

example_SOURCES = $(srcdir)/Example.C \
	$(srcdir)/DennisSchnabel.C $(srcdir)/DennisSchnabel.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H

example_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a

example_CXXFLAGS = $(TRILINOS_CXXFLAGS) \
	-I$(top_srcdir)/src -I$(top_srcdir)/src-epetra

example_LDADD = -L$(top_builddir)/src -L$(top_builddir)/src-epetra \
                -lnox -lnoxepetra 

