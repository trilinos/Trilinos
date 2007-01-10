## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
#            LOCA: Library of Continuation Algorithms Package
#                 Copyright (2005) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# 
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
# Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
# ************************************************************************
# @HEADER

include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(LOCA_EPETRA_TEST_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(LOCA_EPETRA_TEST_INCLUDES))
else
TEST_LIBS     = $(LOCA_EPETRA_TEST_LIBS)
TEST_INCLUDES = $(LOCA_EPETRA_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(TEST_INCLUDES)

EXEEXT = .exe

HOUSEHOLDERBORDEREDSOLVE = householderBorderedSolve
HOUSEHOLDERTRANSBORDEREDSOLVE = householderTransposeBorderedSolve
TRANSPOSESOLVE = transposeSolve
COMPOSITECONSTRAINT = compositeConstraint
COMPOSITECONSTRAINTMVDX = compositeConstraintMVDX
NCRF = naturalContResidualFills

if BUILD_LOCA_EPETRA_TESTS
noinst_PROGRAMS = \
	$(HOUSEHOLDERBORDEREDSOLVE) \
	$(HOUSEHOLDERTRANSBORDEREDSOLVE) \
	$(TRANSPOSESOLVE) \
	$(COMPOSITECONSTRAINT) \
	$(COMPOSITECONSTRAINTMVDX) \
	$(NCRF)
endif

householderBorderedSolve_SOURCES = \
	$(srcdir)/HouseholderBorderedSolve.C 

householderBorderedSolve_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
	$(top_builddir)/src-epetra/libnoxepetra.a \
	$(top_builddir)/test/utils/libnoxtestutils.a \
	$(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.a \
	$(top_builddir)/test/epetra/LOCA_TestProblems/liblocatestproblems.a

householderBorderedSolve_LDADD = $(TEST_LIBS)

householderTransposeBorderedSolve_SOURCES = \
	$(srcdir)/HouseholderTransposeBorderedSolve.C 

householderTransposeBorderedSolve_DEPENDENCIES = \
	$(householderBorderedSolve_DEPENDENCIES)

householderTransposeBorderedSolve_LDADD = \
	$(householderBorderedSolve_LDADD)

transposeSolve_SOURCES = \
	$(srcdir)/TransposeSolve.C 

transposeSolve_DEPENDENCIES = \
	$(householderBorderedSolve_DEPENDENCIES)

transposeSolve_LDADD = \
	$(householderBorderedSolve_LDADD)

compositeConstraint_SOURCES = \
	$(srcdir)/CompositeConstraint.C

compositeConstraint_DEPENDENCIES = \
	$(householderBorderedSolve_DEPENDENCIES)

compositeConstraint_LDADD = \
	$(householderBorderedSolve_LDADD)

compositeConstraintMVDX_SOURCES = \
	$(srcdir)/CompositeConstraintMVDX.C

compositeConstraintMVDX_DEPENDENCIES = \
	$(householderBorderedSolve_DEPENDENCIES)

compositeConstraintMVDX_LDADD = \
	$(householderBorderedSolve_LDADD)

naturalContResidualFills_SOURCES = \
	$(srcdir)/NaturalContResidualFills.C

naturalContResidualFills_DEPENDENCIES = \
	$(householderBorderedSolve_DEPENDENCIES)

naturalContResidualFills_LDADD = \
	$(householderBorderedSolve_LDADD)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)

