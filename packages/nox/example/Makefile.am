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

if BUILD_EXAMPLES

EXAMPLES_LIBRARY = libnoxexample.a
EXAMPLES = rosenbrock broyden

endif

noinst_LIBRARIES = $(EXAMPLES_LIBRARY)
noinst_PROGRAMS = $(EXAMPLES)

# Header files are included here so that they are included in the
# automatically generated distribution.
libnoxexample_a_SOURCES = \
	$(srcdir)/NOX_Example_Vector.H \
	$(srcdir)/NOX_Example_Vector.C \
	$(srcdir)/NOX_Example_Matrix.H \
	$(srcdir)/NOX_Example_Matrix.C \
	$(srcdir)/NOX_Example_Interface.H \
	$(srcdir)/NOX_Example_Group.H \
	$(srcdir)/NOX_Example_Group.C 

libnoxexample_a_CXXFLAGS = -I$(top_srcdir)/src

libnoxexample_a_LIBADD = $(CXXREPOS)

rosenbrock_SOURCES = $(srcdir)/Rosenbrock.C

rosenbrock_DEPENDENCIES = $(top_builddir)/example/libnoxexample.a $(top_builddir)/src/libnox.a 

rosenbrock_CXXFLAGS = -I$(top_srcdir)/src

# Be sure to specify the **exact** library you want since others may be elsewhere in the path
rosenbrock_LDADD = $(top_builddir)/example/libnoxexample.a $(top_builddir)/src/libnox.a 

broyden_SOURCES = $(srcdir)/Broyden.C

broyden_DEPENDENCIES = $(top_builddir)/example/libnoxexample.a $(top_builddir)/src/libnox.a 

broyden_CXXFLAGS = -I$(top_srcdir)/src

# Order matters here!
broyden_LDADD = $(CXXREPOS) -L$(top_builddir)/src -L$(top_builddir)/example -lnoxexample -lnox 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libnoxexample_a_AR = $(ALTERNATE_AR)

else

libnoxexample_a_AR = $(AR) cru

endif

MAINTAINERCLEANFILES = Makefile.in 

MOSTLYCLEANFILES = $(CXXREPOS)


