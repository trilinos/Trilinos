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
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside NOX. The first is
:q# the list of all source and any other header files.
# ------------------------------------------------------------------------

NOX_BLAS_H = \
	$(srcdir)/blas.h \
	$(srcdir)/lapack.h \
	$(srcdir)/NOX_BLAS.H \
	$(srcdir)/NOX_BLAS_Vector.H \
	$(srcdir)/NOX_BLAS_Matrix.H \
	$(srcdir)/NOX_BLAS_Group.H \
	$(srcdir)/NOX_BLAS_Interface.H 

NOX_BLAS = \
	$(srcdir)/NOX_BLAS_Vector.C \
	$(srcdir)/NOX_BLAS_Matrix.C \
	$(srcdir)/NOX_BLAS_Group.C 

LOCA_BLAS_H = \
	$(srcdir)/LOCA_BLAS.H \
	$(srcdir)/LOCA_BLAS_Group.H \
	$(srcdir)/LOCA_BLAS_Interface.H \
	$(srcdir)/LOCA_BLAS_DataOutput.H \
	$(srcdir)/LOCA_BLAS_TPDataOutput.H

LOCA_BLAS = \
	$(srcdir)/LOCA_BLAS_Group.C \
	$(srcdir)/LOCA_BLAS_DataOutput.C \
	$(srcdir)/LOCA_BLAS_TPDataOutput.C 

# ------------------------------------------------------------------------
# Noxblas library specifications
# ------------------------------------------------------------------------

libnoxblas_a_SOURCES = $(NOX_BLAS) $(NOX_BLAS_H)

libnoxblas_a_CXXFLAGS = -I$(top_srcdir)/src-loca/src \
                           -I$(top_srcdir)/src \
                           $(TRILINOS_CXXFLAGS) 

liblocablas_a_SOURCES = $(LOCA_BLAS) $(LOCA_BLAS_H)

liblocablas_a_CXXFLAGS = -I$(top_srcdir)/src-loca/src \
                           -I$(top_srcdir)/src \
                           $(TRILINOS_CXXFLAGS) 

# ------------------------------------------------------------------------
# When to turn it on
# ------------------------------------------------------------------------

if BUILD_LIBLOCABLAS

# List of all libraries to install in $(libexecdir)

lib_LIBRARIES = libnoxblas.a liblocablas.a

# List of all headers to install in $(pkgincludedir)

pkginclude_HEADERS = $(NOX_BLAS_H) $(LOCA_BLAS_H)

endif

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libnoxblas_a_AR = $(ALTERNATE_AR)
liblocablas_a_AR = $(ALTERNATE_AR)

else

libnoxblas_a_AR = $(AR) cru
liblocablas_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libnoxblas_a_LIBADD = $(CXXREPOS)
liblocablas_a_LIBADD = $(CXXREPOS)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






