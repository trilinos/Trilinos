## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
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
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov), Sandia National Laboratories.
# 
# ************************************************************************
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside NOX. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

NOX_TSF_H = \
	$(srcdir)/NOX_TSF.H \
	$(srcdir)/NOX_TSF_Vector.H \
	$(srcdir)/NOX_TSF_Group.H 

NOX_TSF = \
	$(srcdir)/NOX_TSF_Vector.C \
	$(srcdir)/NOX_TSF_Group.C 

# ------------------------------------------------------------------------
# Noxblas library specifications
# ------------------------------------------------------------------------

libnoxtsf_a_SOURCES = $(NOX_TSF) $(NOX_TSF_H) 

libnoxtsf_a_CXXFLAGS = \
	-I$(top_srcdir)/src-loca/src \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/../TSFExtended/src/LinearSolvers \
	-I$(top_srcdir)/../TSFExtended/src/Preconditioners \
	-I$(top_srcdir)/../TSFExtended/src/Operators \
	-I$(top_srcdir)/../TSFExtended/src/Utils \
	-I$(top_srcdir)/../TSFExtended/src/Vectors \
	-I$(top_srcdir)/../TSFExtended/src/Concrete \
	-I$(top_srcdir)/../TSFCoreUtils/src \
	-I$(top_srcdir)/../TSFCoreUtils/src/RTOpStdOps \
	-I$(top_srcdir)/../TSFCoreUtils/src/RTOpPack \
	-I$(top_srcdir)/../TSFCoreUtils/src/MoochoUtilities \
	-I$(top_srcdir)/../teuchos/src \
	-I$(top_builddir)/../teuchos/src \
	-I$(top_builddir)/../TSFExtended/src \
	-I$(top_builddir)/../TSFCoreUtils/src \
	-I$(top_builddir)/../TSFCore/src \
	-I$(top_srcdir)/../TSFCore \
	-I$(top_srcdir)/../TSFCore/src \
	-I$(top_srcdir)/../TSFCore/src/adapters/epetra \
	-I$(top_srcdir)/../TSFCore/src/adapters/mpi-base \
	-I$(top_srcdir)/../TSFCore/src/interfaces/Core \
	$(TRILINOS_CXXFLAGS) 

# ------------------------------------------------------------------------
# When to turn it on
# ------------------------------------------------------------------------

if BUILD_LIBNOXTSF

# List of all libraries to install in $(libexecdir)

lib_LIBRARIES = libnoxtsf.a

# List of all headers to install in $(includedir)

include_HEADERS = $(NOX_TSF_H) 

endif

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libnoxtsf_a_AR = $(ALTERNATE_AR)

else

libnoxtsf_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libnoxtsf_a_LIBADD = $(CXXREPOS)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






