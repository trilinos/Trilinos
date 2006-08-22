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


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside NOX. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

LOCA_TESTPROBLEMS_H = \
	$(srcdir)/Basis.H \
	$(srcdir)/FiniteElementProblem.H \
	$(srcdir)/Problem_Interface.H \
	$(srcdir)/Tcubed_FiniteElementProblem.H \
	$(srcdir)/Pitchfork_FiniteElementProblem.H \
	$(srcdir)/NormConstraint.H \
	$(srcdir)/LinearConstraint.H 

LOCA_TESTPROBLEMS = \
	$(srcdir)/Basis.C \
	$(srcdir)/Problem_Interface.C \
	$(srcdir)/Tcubed_FiniteElementProblem.C \
	$(srcdir)/Pitchfork_FiniteElementProblem.C \
	$(srcdir)/NormConstraint.C \
	$(srcdir)/LinearConstraint.C

# ------------------------------------------------------------------------
# Include NOX export makefile
# ------------------------------------------------------------------------
include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(LOCA_EPETRA_TEST_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(LOCA_EPETRA_TEST_INCLUDES))
else
TEST_LIBS     = $(LOCA_EPETRA_TEST_LIBS)
TEST_INCLUDES = $(LOCA_EPETRA_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(TEST_INCLUDES)

# ------------------------------------------------------------------------
# locatestproblems library specifications
# ------------------------------------------------------------------------

liblocatestproblems_a_SOURCES = $(LOCA_TESTPROBLEMS) $(LOCA_TESTPROBLEMS_H)

if BUILD_LOCA_EPETRA_TESTS
noinst_LIBRARIES = liblocatestproblems.a
endif

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

liblocatestproblems_a_AR = $(ALTERNATE_AR)

else

liblocatestproblems_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

liblocatestproblems_a_LIBADD = $(CXXREPOS) $(XTRALDADD)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)



