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

if BUILD_LOCA_LAPACK_EXAMPLES

PITCHFORK = pitchfork
PITCHFORKCONT = pitchforkContinuation

endif

noinst_PROGRAMS = $(PITCHFORK) $(PITCHFORKCONT)

pitchfork_SOURCES = $(srcdir)/Pitchfork.C \
		$(srcdir)/PitchforkProblemInterface.H \
		$(srcdir)/PitchforkProblemInterface.C

pitchforkContinuation_SOURCES = $(srcdir)/PitchforkContinuation.C \
		$(srcdir)/PitchforkProblemInterface.H \
		$(srcdir)/PitchforkProblemInterface.C

pitchfork_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
	$(top_builddir)/src-lapack/libnoxlapack.a \
	$(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-lapack/liblocalapack.a

pitchforkContinuation_DEPENDENCIES = $(pitchfork_DEPENDENCIES)

pitchfork_CXXFLAGS = $(TRILINOS_CXXFLAGS) \
		   -I$(top_srcdir)/src \
		   -I$(top_srcdir)/src-lapack \
		   -I$(top_srcdir)/src-loca/src \
                   -I$(top_srcdir)/src-loca/src-lapack

pitchforkContinuation_CXXFLAGS = $(pitchfork_CXXFLAGS)

pitchfork_LDADD = \
		-L$(top_builddir)/src \
		-L$(top_builddir)/src-lapack \
                -L$(top_builddir)/src-loca/src \
                -L$(top_builddir)/src-loca/src-lapack \
                -llocalapack -lloca -lnoxlapack -lnox

pitchforkContinuation_LDADD = $(pitchfork_LDADD)
