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

noinst_PROGRAMS = rosenbrock broyden

rosenbrock_SOURCES = \
                   $(srcdir)/Rosenbrock.C 

rosenbrock_DEPENDENCIES = $(top_builddir)/src/libnox.a \
                          $(top_builddir)/src-lapack/libnoxlapack.a \
                          $(top_builddir)/test/utils/libnoxtestutils.a

rosenbrock_CXXFLAGS = -I$(top_srcdir)/src \
                      -I$(top_srcdir)/src-lapack \
                      -I$(top_srcdir)/test/utils \
                      $(TRILINOS_CXXFLAGS)

rosenbrock_LDADD = $(top_builddir)/src-lapack/libnoxlapack.a \
		   $(top_builddir)/src/libnox.a \
		   $(top_builddir)/test/utils/libnoxtestutils.a \
		   $(LIBS) 

broyden_SOURCES = $(srcdir)/Broyden.C

broyden_DEPENDENCIES = $(rosenbrock_DEPENDENCIES)

broyden_CXXFLAGS = $(rosenbrock_CXXFLAGS)

broyden_LDADD = $(rosenbrock_LDADD)
