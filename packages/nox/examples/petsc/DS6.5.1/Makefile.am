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
#@HEADER

if BUILD_PETSC_EXAMPLES
include ${PETSC_DIR}/bmake/common/base
endif

EXEEXT = .exe

if BUILD_PETSC_EXAMPLES
noinst_PROGRAMS = example
endif

example_SOURCES = $(srcdir)/Example.C \
	$(srcdir)/DennisSchnabel.C $(srcdir)/DennisSchnabel.H \
	$(srcdir)/Problem_Interface.C $(srcdir)/Problem_Interface.H

example_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-petsc/libnoxpetsc.a

example_CXXFLAGS = $(PETSC_CXXFLAGS) \
	-I$(top_srcdir)/src -I$(top_srcdir)/src-petsc

example_LDADD = -L$(top_builddir)/src -L$(top_builddir)/src-petsc \
                -lnoxpetsc ${PETSC_LIB} -lnox

