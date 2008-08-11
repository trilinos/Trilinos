#@HEADER
#  ************************************************************************
#
#                Kokkos: Fast Core Kernels Package
#               Copyright (2004) Sandia Corporation
#
#  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
#  license for use of this work by or on behalf of the U.S. Government.
#
#  This library is free software; you can redistribute it and or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation; either version 2.1 of the
#  License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
#  USA
#  Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
#  ************************************************************************
#@HEADER


EXEEXT = .exe

noinst_PROGRAMS = SparseMatByMultiVec5_example 

SparseMatByMultiVec5_example_SOURCES = $(srcdir)/csrcsc.f $(srcdir)/scscmm1a.f \
     $(srcdir)/scscmm2a.f $(srcdir)/scscrd.f $(srcdir)/scsrmm1a5.f \
     $(srcdir)/scsrmm2a5.f $(srcdir)/snrm2.f $(srcdir)/srelerr.f \
     $(srcdir)/swaxpy.f $(srcdir)/timer.c $(srcdir)/tsploops.f

#This is a stand-alone example
#SparseMatByMultiVec5_example_DEPENDENCIES = $(top_builddir)/src/libkokkos.a

#SparseMatByMultiVec5_example_CXXFLAGS = -I$(top_srcdir)/src

#SparseMatByMultiVec5_example_LDADD = $(top_builddir)/src/libkokkos.a \
#	$(LIBS) $(BLAS_LIBS) $(FLIBS)

