## $Id$ 
## $Source$ 
#@HEADER
# ************************************************************************
# 
#                  LOCA Continuation Algorithm Package
#                 Copyright (2005) Sandia Corporation
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

if BUILD_LOCA_LAPACK_EXAMPLES

BRUSSELATORCONT = brusselatorContinuation
BRUSSELATORHOPF = brusselatorHopf
BRUSSELATORHOPFCONTINUATION = brusselatorHopfContinuation

endif

noinst_PROGRAMS = $(BRUSSELATORCONT) $(BRUSSELATORHOPF) $(BRUSSELATORHOPFCONTINUATION)

brusselatorContinuation_SOURCES = $(srcdir)/BrusselatorContinuation.C \
		$(srcdir)/BrusselatorProblemInterface.H \
		$(srcdir)/BrusselatorProblemInterface.C

brusselatorHopf_SOURCES = $(srcdir)/BrusselatorHopf.C \
		$(srcdir)/BrusselatorProblemInterface.H \
		$(srcdir)/BrusselatorProblemInterface.C

brusselatorHopfContinuation_SOURCES = $(srcdir)/BrusselatorHopfContinuation.C \
		$(srcdir)/BrusselatorProblemInterface.H \
		$(srcdir)/BrusselatorProblemInterface.C

brusselatorContinuation_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
	$(top_builddir)/src-lapack/libnoxlapack.a \
	$(top_builddir)/src-loca/src/libloca.la \
        $(top_builddir)/src-loca/src-lapack/liblocalapack.la

brusselatorHopf_DEPENDENCIES = $(brusselatorContinuation_DEPENDENCIES)

brusselatorHopfContinuation_DEPENDENCIES = $(brusselatorContinuation_DEPENDENCIES)

brusselatorContinuation_CXXFLAGS = $(TRILINOS_CXXFLAGS) \
		   -I$(top_srcdir)/src \
		   -I$(top_srcdir)/src-lapack \
		   -I$(top_srcdir)/src-loca/src \
                   -I$(top_srcdir)/src-loca/src-lapack

brusselatorHopf_CXXFLAGS = $(brusselatorContinuation_CXXFLAGS)

brusselatorHopfContinuation_CXXFLAGS = $(brusselatorContinuation_CXXFLAGS)

brusselatorContinuation_LDADD = \
	$(top_builddir)/src-loca/src-lapack/liblocalapack.la \
	$(top_builddir)/src-loca/src/libloca.la \
	$(top_builddir)/src-lapack/libnoxlapack.a \
	$(top_builddir)/src/libnox.a \
	$(LIBS)

brusselatorHopf_LDADD = $(brusselatorContinuation_LDADD)

brusselatorHopfContinuation_LDADD = $(brusselatorContinuation_LDADD)
