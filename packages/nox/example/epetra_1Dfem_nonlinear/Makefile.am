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

example_CXXFLAGS = -DEPETRA_$(TRILINOS_COMM) \
	-I$(TRILINOS_HOME)/packages/epetra/src \
	-I$(top_srcdir)/src

example_LDADD = \
	-L$(top_builddir)/src \
	-lnox -lnoxepetra \
	-L$(TRILINOS_HOME)/lib/$(TRILINOS_TARGET) \
	-lepetra -laztecoo -lifpack \
	-L/usr/lib -L/usr/lib/local \
	-llapack -lspblas -lblas -ly12m -lg2c -lm

