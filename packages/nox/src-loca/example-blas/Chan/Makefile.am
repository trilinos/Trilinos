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
#@HEADER

if BUILD_LOCA_BLAS_EXAMPLES

CHANCONT = chanContinuation
CHANTP = chanTP
CHANTPCONT = chanTPContinuation

endif

noinst_PROGRAMS = $(CHANCONT) $(CHANTP) $(CHANTPCONT)

chanContinuation_SOURCES = $(srcdir)/ChanContinuation.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTP_SOURCES = $(srcdir)/ChanTP.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTPContinuation_SOURCES = $(srcdir)/ChanTPContinuation.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanContinuation_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
	$(top_builddir)/src-loca/src/libloca.a \
	$(top_builddir)/src-loca/src-blas/libnoxblas.a \
        $(top_builddir)/src-loca/src-blas/liblocablas.a

chanTP_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanTPContinuation_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanContinuation_CXXFLAGS = $(TRILINOS_CXXFLAGS) \
		   -I$(top_srcdir)/src \
		   -I$(top_srcdir)/src-loca/src \
                   -I$(top_srcdir)/src-loca/src-blas

chanTP_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanTPContinuation_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanContinuation_LDADD = \
		-L$(top_builddir)/src \
                -L$(top_builddir)/src-loca/src \
                -L$(top_builddir)/src-loca/src-blas \
                -llocablas -lloca -lnoxblas -lnox

chanTP_LDADD = $(chanContinuation_LDADD)

chanTPContinuation_LDADD = $(chanContinuation_LDADD)