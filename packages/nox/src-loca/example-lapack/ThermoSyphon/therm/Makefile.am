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
# therm library specifications
# ------------------------------------------------------------------------

libtherm_a_SOURCES = \
		$(srcdir)/array.h \
		$(srcdir)/limits.h \
		$(srcdir)/indata.h \
		$(srcdir)/precond.h \
		$(srcdir)/3579.h \
		$(srcdir)/indata.f \
		$(srcdir)/full_pack.f \
		$(srcdir)/compB.f \
		$(srcdir)/compR.f \
		$(srcdir)/compI.f \
		$(srcdir)/convol.f \
		$(srcdir)/icostr.f \
		$(srcdir)/initcos.f \
		$(srcdir)/findlargest.f \
		$(srcdir)/therminit.f \
		$(srcdir)/thermfill.f

THERM = therm

# ------------------------------------------------------------------------
# When to turn it on
# ------------------------------------------------------------------------

if BUILD_LOCA_LAPACK_EXAMPLES

# List of all libraries to install in $(libexecdir)

lib_LIBRARIES = libtherm.a

endif

noinst_PROGRAMS = $(THERM)

therm_SOURCES = $(srcdir)/main.f

therm_DEPENDENCIES = \
	$(top_builddir)/src-loca/example-lapack/ThermoSyphon/therm/libtherm.a

therm_LDADD = \
	-L$(top_builddir)/src-loca/example-lapack/ThermoSyphon/therm \
	-ltherm
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libtherm_a_AR = $(ALTERNATE_AR)

else

libtherm_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libtherm_a_LIBADD = $(CXXREPOS)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

