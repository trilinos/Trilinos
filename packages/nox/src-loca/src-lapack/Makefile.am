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

LOCA_LAPACK_H = \
	$(srcdir)/LOCA_LAPACK.H \
	$(srcdir)/LOCA_LAPACK_Group.H \
	$(srcdir)/LOCA_LAPACK_Interface.H \
	$(srcdir)/LOCA_Eigensolver_DGGEVStrategy.H \
	$(srcdir)/LOCA_BorderedSystem_LAPACKDirectSolve.H \
	$(srcdir)/LOCA_LAPACK_Factory.H

LOCA_LAPACK = \
	$(srcdir)/LOCA_LAPACK_Group.C \
	$(srcdir)/LOCA_Eigensolver_DGGEVStrategy.C \
	$(srcdir)/LOCA_BorderedSystem_LAPACKDirectSolve.C \
	$(srcdir)/LOCA_LAPACK_Factory.C

# ------------------------------------------------------------------------
# localapack library specifications
# ------------------------------------------------------------------------

liblocalapack_la_SOURCES = $(LOCA_LAPACK) $(LOCA_LAPACK_H)

liblocalapack_la_CXXFLAGS = \
	-I$(top_srcdir)/src-loca/src \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/src-lapack \
	$(TRILINOS_CXXFLAGS) $(LIBTOOL_COMPILER_FLAGS)

liblocalapack_la_LIBADD = \
	$(top_builddir)/src-loca/src/libloca.la \
	$(CXXREPOS) $(XTRALDADD) $(LIBTOOL_CCLINKER_FLAGS)

# ------------------------------------------------------------------------
# When to turn it on
# ------------------------------------------------------------------------

if BUILD_LIBLOCALAPACK

# List of all libraries to install in $(libexecdir)

lib_LTLIBRARIES = liblocalapack.la

# List of all headers to install in $(includedir)

include_HEADERS = $(LOCA_LAPACK_H)

endif

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





