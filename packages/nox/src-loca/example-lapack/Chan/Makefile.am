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

if BUILD_LOCA_LAPACK_EXAMPLES

CHANCONT = chanContinuation
CHANTP = chanTP
CHANTPCONT = chanTPContinuation
CHANHOMOTOPY = chanHomotopy
CHANTPMOD = chanTP_ModifiedBordering
CHANTPCONTMOD = chanTPContinuation_ModifiedBordering
CHANTPNDMOD = chanTP_NicDayModifiedBordering
CHANTPCONTNDMOD = chanTPContinuation_NicDayModifiedBordering

if BUILD_LIBLOCAMF
CHANMULTICONT = chanMultiContinuation
endif

endif

noinst_PROGRAMS = $(CHANCONT) $(CHANTP) $(CHANTPCONT) $(CHANHOMOTOPY) \
	          $(CHANTPMOD) $(CHANTPCONTMOD) $(CHANTPNDMOD) \
		  $(CHANTPCONTNDMOD) $(CHANMULTICONT)

chanContinuation_SOURCES = $(srcdir)/ChanContinuation.C \
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
                   -I$(top_srcdir)/src-loca/src-lapack

chanContinuation_LDADD = \
		-L$(top_builddir)/src \
		-L$(top_builddir)/src-lapack \
                -L$(top_builddir)/src-loca/src \
                -L$(top_builddir)/src-loca/src-lapack \
                -llocalapack -lloca -lnoxlapack -lnox $(LIBS)

chanTP_SOURCES = $(srcdir)/ChanTP.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTP_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanTP_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanTP_LDADD = $(chanContinuation_LDADD)

chanTPContinuation_SOURCES = $(srcdir)/ChanTPContinuation.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTPContinuation_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanTPContinuation_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanTPContinuation_LDADD = $(chanContinuation_LDADD)

chanHomotopy_SOURCES = $(srcdir)/ChanHomotopy.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanHomotopy_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanHomotopy_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanHomotopy_LDADD = $(chanContinuation_LDADD)

chanTP_ModifiedBordering_SOURCES = $(srcdir)/ChanTP_ModifiedBordering.C \
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

chanTP_NicDayModifiedBordering_SOURCES = \
		$(srcdir)/ChanTP_NicDayModifiedBordering.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTP_NicDayModifiedBordering_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanTP_NicDayModifiedBordering_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanTP_NicDayModifiedBordering_LDADD = $(chanContinuation_LDADD)

chanTPContinuation_NicDayModifiedBordering_SOURCES = \
		$(srcdir)/ChanTPContinuation_NicDayModifiedBordering.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTPContinuation_NicDayModifiedBordering_DEPENDENCIES = \
	$(chanContinuation_DEPENDENCIES)

chanTPContinuation_NicDayModifiedBordering_CXXFLAGS = $(chanContinuation_CXXFLAGS)

chanTPContinuation_NicDayModifiedBordering_LDADD = $(chanContinuation_LDADD)

chanMultiContinuation_SOURCES = $(srcdir)/ChanMultiContinuation.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanMultiContinuation_DEPENDENCIES = $(chanContinuation_DEPENDENCIES) \
                                 $(top_builddir)/src-loca/src-mf/liblocamf.a \
                                 $(MF_LIBDIR)/libExpCmp.a \
                                 $(MF_LIBDIR)/libMF.a

chanMultiContinuation_CXXFLAGS = $(chanContinuation_CXXFLAGS) \
                                    -I$(top_srcdir)/src-loca/src-mf \
                                    -I$(MF_INCDIR)

chanMultiContinuation_LDADD = -L$(top_builddir)/src-loca/src-mf -llocamf \
                              -L$(MF_LIBDIR) -lMF -lExpCmp \
                                 $(chanContinuation_LDADD)
