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

#EXEEXT = .exe

CHANCONT = chanContinuation_exe
CHANTP = chanTP_exe
CHANTPCONT = chanTPContinuation_exe
CHANHOMOTOPY = chanHomotopy_exe
CHANTPMOD = chanTP_ModifiedBordering_exe
CHANTPCONTMOD = chanTPContinuation_ModifiedBordering_exe

noinst_PROGRAMS = $(CHANCONT) $(CHANTP) $(CHANTPCONT) $(CHANHOMOTOPY) \
	          $(CHANTPMOD) $(CHANTPCONTMOD)

chanContinuation_exe_SOURCES = \
                $(srcdir)/ChanContinuation.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanContinuation_exe_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
	$(top_builddir)/src-lapack/libnoxlapack.a \
	$(top_builddir)/src-loca/src/libloca.la \
        $(top_builddir)/src-loca/src-lapack/liblocalapack.la

chanContinuation_exe_CXXFLAGS = $(TRILINOS_CXXFLAGS) \
		   -I$(top_srcdir)/src \
		   -I$(top_srcdir)/src-lapack \
		   -I$(top_srcdir)/src-loca/src \
                   -I$(top_srcdir)/src-loca/src-lapack \
                   -I$(top_srcdir)/test/utils \
	           $(TRILINOS_CXXFLAGS)

chanContinuation_exe_LDADD = \
                $(top_builddir)/src-loca/src-lapack/liblocalapack.la \
                $(top_builddir)/src-loca/src/libloca.la \
		-L$(top_builddir)/src \
		-L$(top_builddir)/src-lapack \
                -lnoxlapack -lnox \
                $(LIBS) $(FLIBS)

chanTP_exe_SOURCES = \
                $(srcdir)/ChanTP.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTP_exe_DEPENDENCIES = $(chanContinuation_exe_DEPENDENCIES)

chanTP_exe_CXXFLAGS = $(chanContinuation_exe_CXXFLAGS)

chanTP_exe_LDADD = $(chanContinuation_exe_LDADD)

chanTPContinuation_exe_SOURCES = \
                $(srcdir)/ChanTPContinuation.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTPContinuation_exe_DEPENDENCIES = $(chanContinuation_exe_DEPENDENCIES)

chanTPContinuation_exe_CXXFLAGS = $(chanContinuation_exe_CXXFLAGS)

chanTPContinuation_exe_LDADD = $(chanContinuation_exe_LDADD)

chanHomotopy_exe_SOURCES = \
                $(srcdir)/ChanHomotopy.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanHomotopy_exe_DEPENDENCIES = $(chanContinuation_exe_DEPENDENCIES)

chanHomotopy_exe_CXXFLAGS = $(chanContinuation_exe_CXXFLAGS)

chanHomotopy_exe_LDADD = $(chanContinuation_exe_LDADD)

chanTP_ModifiedBordering_exe_SOURCES = \
                $(srcdir)/ChanTP_ModifiedBordering.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTP_ModifiedBordering_exe_DEPENDENCIES = $(chanContinuation_exe_DEPENDENCIES)

chanTP_ModifiedBordering_exe_CXXFLAGS = $(chanContinuation_exe_CXXFLAGS)

chanTP_ModifiedBordering_exe_LDADD = $(chanContinuation_exe_LDADD)

chanTPContinuation_ModifiedBordering_exe_SOURCES = \
		$(srcdir)/ChanTPContinuation_ModifiedBordering.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTPContinuation_ModifiedBordering_exe_DEPENDENCIES = \
	$(chanContinuation_exe_DEPENDENCIES)

chanTPContinuation_ModifiedBordering_exe_CXXFLAGS = $(chanContinuation_exe_CXXFLAGS)

chanTPContinuation_ModifiedBordering_exe_LDADD = $(chanContinuation_exe_LDADD)
