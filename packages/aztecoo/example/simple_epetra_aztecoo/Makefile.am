#@HEADER
# ************************************************************************
# 
#        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = epetra_crsmatrix

epetra_crsmatrix_SOURCES = $(srcdir)/cxx_main.cpp

#epetra_crsmatrix_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a $(top_builddir)/../epetra/src/libepetra.a

epetra_crsmatrix_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/../epetra/src \
	-I$(top_builddir)/../epetra/src

epetra_crsmatrix_LDADD = $(top_builddir)/src/libaztecoo.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)


