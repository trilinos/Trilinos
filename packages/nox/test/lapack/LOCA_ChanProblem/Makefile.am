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

EXEEXT = .exe

CHANCONT = chanContinuation
CHANTP = chanTP
CHANTPCONT = chanTPContinuation
CHANHOMOTOPY = chanHomotopy
CHANTPMOD = chanTP_ModifiedBordering
CHANTPCONTMOD = chanTPContinuation_ModifiedBordering

noinst_PROGRAMS = $(CHANCONT) $(CHANTP) $(CHANTPCONT) $(CHANHOMOTOPY) \
	          $(CHANTPMOD) $(CHANTPCONTMOD)

chanContinuation_SOURCES = \
                $(srcdir)/ChanContinuation.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanContinuation_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
	$(top_builddir)/src-lapack/libnoxlapack.a \
	$(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-lapack/liblocalapack.a

chanContinuation_CXXFLAGS = $(TRILINOS_CXXFLAGS) \
		   -I$(top_srcdir)/src \
		   -I$(top_srcdir)/src-lapack \
		   -I$(top_srcdir)/src-loca/src \
                   -I$(top_srcdir)/src-loca/src-lapack \
                   -I$(top_srcdir)/test/utils

chanContinuation_LDADD = \
		-L$(top_builddir)/src \
		-L$(top_builddir)/src-lapack \
                -L$(top_builddir)/src-loca/src \
                -L$(top_builddir)/src-loca/src-lapack \
                -llocalapack -lloca -lnoxlapack -lnox

chanTP_SOURCES = \
                $(srcdir)/ChanTP.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTP_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanTP_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanTP_LDADD = $(chanContinuation_LDADD)

chanTPContinuation_SOURCES = \
                $(srcdir)/ChanTPContinuation.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTPContinuation_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanTPContinuation_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanTPContinuation_LDADD = $(chanContinuation_LDADD)

chanHomotopy_SOURCES = \
                $(srcdir)/ChanHomotopy.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanHomotopy_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanHomotopy_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanHomotopy_LDADD = $(chanContinuation_LDADD)

chanTP_ModifiedBordering_SOURCES = \
                $(srcdir)/ChanTP_ModifiedBordering.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTP_ModifiedBordering_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanTP_ModifiedBordering_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanTP_ModifiedBordering_LDADD = $(chanContinuation_LDADD)

chanTPContinuation_ModifiedBordering_SOURCES = \
		$(srcdir)/ChanTPContinuation_ModifiedBordering.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTPContinuation_ModifiedBordering_DEPENDENCIES = \
	$(chanContinuation_DEPENDENCIES)

chanTPContinuation_ModifiedBordering_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanTPContinuation_ModifiedBordering_LDADD = $(chanContinuation_LDADD)