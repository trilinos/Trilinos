#@HEADER
# ************************************************************************
# 
#          Tpetra: Templated Linear Algebra Services Package
#                 Copyright (2004) Sandia Corporation
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = CisMatrix_test

CisMatrix_test_SOURCES = \
	$(srcdir)/cxx_main.cpp

CisMatrix_test_DEPENDENCIES = \
	$(top_builddir)/src/libtpetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../kokkos/src/libkokkos.a 

CisMatrix_test_CXXFLAGS = \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/../teuchos/src \
	-I$(top_srcdir)/../kokkos/src \
	-I$(top_builddir)/../teuchos/src \
	-I$(top_builddir)/../kokkos/src

CisMatrix_test_LDADD = \
	$(top_builddir)/src/libtpetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../kokkos/src/libkokkos.a \
	$(LIBS) $(BLAS_LIBS) $(FLIBS)

