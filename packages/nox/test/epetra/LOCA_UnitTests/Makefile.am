## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                  LOCA Continuation Algorithm Package
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
# Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
# (etphipp@sandia.gov), Sandia National Laboratories.
# 
# ************************************************************************
# @HEADER

HOUSEHOLDERBORDEREDSOLVE = householderBorderedSolve_exe
HOUSEHOLDERBORDEREDSOLVECONTIGUOUS = householderBorderedSolveContiguous_exe
COMPOSITECONSTRAINT = compositeConstraint_exe
COMPOSITECONSTRAINTMVDX = compositeConstraintMVDX_exe

noinst_PROGRAMS = \
	$(HOUSEHOLDERBORDEREDSOLVE) \
	$(HOUSEHOLDERBORDEREDSOLVECONTIGUOUS) \
	$(COMPOSITECONSTRAINT) \
	$(COMPOSITECONSTRAINTMVDX)

householderBorderedSolve_exe_SOURCES = \
	$(srcdir)/HouseholderBorderedSolve.C 

householderBorderedSolve_exe_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
	$(top_builddir)/src-epetra/libnoxepetra.a \
	$(top_builddir)/test/utils/libnoxtestutils.a \
	$(top_builddir)/src-loca/src/libloca.la \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.la \
	$(top_builddir)/test/epetra/LOCA_TestProblems/liblocatestproblems.la

householderBorderedSolve_exe_CXXFLAGS = \
	$(TRILINOS_CXXFLAGS) \
	$(LIBTOOL_COMPILER_FLAGS) \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/src-epetra \
	-I$(top_srcdir)/src-loca/src \
        -I$(top_srcdir)/src-loca/src-epetra \
        -I$(top_srcdir)/test/utils \
	-I$(top_srcdir)/test/epetra/LOCA_TestProblems

householderBorderedSolve_exe_LDADD = \
	$(top_builddir)/test/epetra/LOCA_TestProblems/liblocatestproblems.la \
        $(top_builddir)/src-loca/src-epetra/liblocaepetra.la \
        $(top_builddir)/src-loca/src/libloca.la \
	-L$(top_builddir)/src \
	-L$(top_builddir)/src-epetra \
	-L$(top_builddir)/test/utils \
        -lnoxtestutils -lnoxepetra -lnox \
        $(XTRALDADD) $(LIBTOOL_CCLINKER_FLAGS) $(LIBS)

householderBorderedSolveContiguous_exe_SOURCES = \
	$(srcdir)/HouseholderBorderedSolveContiguous.C 

householderBorderedSolveContiguous_exe_DEPENDENCIES = \
	$(householderBorderedSolve_exe_DEPENDENCIES)

householderBorderedSolveContiguous_exe_CXXFLAGS = \
	$(householderBorderedSolve_exe_CXXFLAGS)

householderBorderedSolveContiguous_exe_LDADD = \
	$(householderBorderedSolve_exe_LDADD)

compositeConstraint_exe_SOURCES = \
	$(srcdir)/CompositeConstraint.C

compositeConstraint_exe_DEPENDENCIES = \
	$(householderBorderedSolve_exe_DEPENDENCIES)

compositeConstraint_exe_CXXFLAGS = \
	$(householderBorderedSolve_exe_CXXFLAGS)

compositeConstraint_exe_LDADD = \
	$(householderBorderedSolve_exe_LDADD)

compositeConstraintMVDX_exe_SOURCES = \
	$(srcdir)/CompositeConstraintMVDX.C

compositeConstraintMVDX_exe_DEPENDENCIES = \
	$(householderBorderedSolve_exe_DEPENDENCIES)

compositeConstraintMVDX_exe_CXXFLAGS = \
	$(householderBorderedSolve_exe_CXXFLAGS)

compositeConstraintMVDX_exe_LDADD = \
	$(householderBorderedSolve_exe_LDADD)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)

# ------------------------------------------------------------------------
# Include global libtool definitions
# ------------------------------------------------------------------------
include $(top_srcdir)/config/libtool_macros.am
