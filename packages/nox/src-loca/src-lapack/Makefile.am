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

LOCA_LAPACK_H = \
	$(srcdir)/LOCA_LAPACK.H \
	$(srcdir)/LOCA_LAPACK_Group.H \
	$(srcdir)/LOCA_LAPACK_Interface.H \
	$(srcdir)/LOCA_Eigensolver_DGGEVStrategy.H \
	$(srcdir)/LOCA_BorderedSolver_LAPACKDirectSolve.H \
	$(srcdir)/LOCA_LAPACK_Factory.H

LOCA_LAPACK = \
	$(srcdir)/LOCA_LAPACK_Group.C \
	$(srcdir)/LOCA_Eigensolver_DGGEVStrategy.C \
	$(srcdir)/LOCA_BorderedSolver_LAPACKDirectSolve.C \
	$(srcdir)/LOCA_LAPACK_Factory.C

# ------------------------------------------------------------------------
# localapack library specifications
# ------------------------------------------------------------------------

liblocalapack_a_SOURCES = $(LOCA_LAPACK) $(LOCA_LAPACK_H)

# ------------------------------------------------------------------------
# Includes for required and optional library dependencies
# ------------------------------------------------------------------------
include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(LOCA_INCLUDES))
else
EXPORT_INCLUDES = $(LOCA_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# When to turn it on
# ------------------------------------------------------------------------

if BUILD_LIBLOCALAPACK

# List of all libraries to install in $(libexecdir)

lib_LIBRARIES = liblocalapack.a

# List of all headers to install in $(includedir)

include_HEADERS = $(LOCA_LAPACK_H)

endif

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

liblocalapack_a_AR = $(ALTERNATE_AR)

else

liblocalapack_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

liblocalapack_a_LIBADD = $(CXXREPOS) $(XTRALDADD)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)



