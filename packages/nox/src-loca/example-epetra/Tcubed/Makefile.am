## $Id$ 
## $Source$ 
#@HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#   
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#   
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov).
# 
# ************************************************************************
#@HEADER

if BUILD_LOCA_EPETRA_EXAMPLES

EXAMPLE = example
EXAMPLETP = exampleTP
EXAMPLETPCONT = exampleTPContinuation
EXAMPLETPMOD = exampleTP_ModifiedBordering
EXAMPLETPCONTMOD = exampleTPContinuation_ModifiedBordering

endif

noinst_PROGRAMS = $(EXAMPLE) $(EXAMPLETP) $(EXAMPLETPCONT) $(EXAMPLETPMOD) \
		  $(EXAMPLETPCONTMOD)

example_SOURCES = $(srcdir)/Example.C \
	$(srcdir)/FiniteElementProblem.C $(srcdir)/FiniteElementProblem.H \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H

example_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a \
	$(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a

example_CXXFLAGS = $(TRILINOS_CXXFLAGS) \
		   -I$(top_srcdir)/src \
                   -I$(top_srcdir)/src-epetra \
		   -I$(top_srcdir)/src-loca/src \
                   -I$(top_srcdir)/src-loca/src-epetra

example_LDADD = -L$(top_builddir)/src \
                -L$(top_builddir)/src-epetra \
                -L$(top_builddir)/src-loca/src \
                -L$(top_builddir)/src-loca/src-epetra \
		-llocaepetra -lnoxepetra -lloca -lnox $(LIBS)

exampleTP_SOURCES = $(srcdir)/ExampleTP.C \
	$(srcdir)/FiniteElementProblem.C $(srcdir)/FiniteElementProblem.H \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H

exampleTP_DEPENDENCIES = $(example_DEPENDENCIES)

exampleTP_CXXFLAGS = $(example_CXXFLAGS)

exampleTP_LDADD = $(example_LDADD)

exampleTPContinuation_SOURCES = $(srcdir)/ExampleTPContinuation.C \
	$(srcdir)/FiniteElementProblem.C $(srcdir)/FiniteElementProblem.H \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H

exampleTPContinuation_DEPENDENCIES = $(example_DEPENDENCIES)

exampleTPContinuation_CXXFLAGS = $(example_CXXFLAGS)

exampleTPContinuation_LDADD = $(example_LDADD)

exampleTP_ModifiedBordering_SOURCES = $(srcdir)/ExampleTP_ModifiedBordering.C \
	$(srcdir)/FiniteElementProblem.C $(srcdir)/FiniteElementProblem.H \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H

exampleTP_ModifiedBordering_DEPENDENCIES = $(example_DEPENDENCIES)

exampleTP_ModifiedBordering_CXXFLAGS = $(example_CXXFLAGS)

exampleTP_ModifiedBordering_LDADD = $(example_LDADD)

exampleTPContinuation_ModifiedBordering_SOURCES = \
	$(srcdir)/ExampleTPContinuation_ModifiedBordering.C \
	$(srcdir)/FiniteElementProblem.C $(srcdir)/FiniteElementProblem.H \
	$(srcdir)/Basis.C $(srcdir)/Basis.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H

exampleTPContinuation_ModifiedBordering_DEPENDENCIES = $(example_DEPENDENCIES)

exampleTPContinuation_ModifiedBordering_CXXFLAGS = $(example_CXXFLAGS)

exampleTPContinuation_ModifiedBordering_LDADD = $(example_LDADD)
