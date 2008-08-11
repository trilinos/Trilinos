## $Id$ 
## $Source$ 
#@HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
#            LOCA: Library of Continuation Algorithms Package
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
# 
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
# Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
# ************************************************************************
#@HEADER

# Include NOX export makefile
include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(LOCA_LAPACK_TEST_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(LOCA_LAPACK_TEST_INCLUDES))
else
TEST_LIBS     = $(LOCA_LAPACK_TEST_LIBS)
TEST_INCLUDES = $(LOCA_LAPACK_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(TEST_INCLUDES) -I$(top_srcdir)/test/lapack/LOCA_TestProblems

EXEEXT = .exe

CHANCONT = chanContinuation
CHANTPCONT = chanTPContinuation
CHANHOMOTOPY = chanHomotopy

if BUILD_LOCA_LAPACK_TESTS
noinst_PROGRAMS = $(CHANCONT) $(CHANTPCONT) $(CHANHOMOTOPY)
endif

chanContinuation_SOURCES = $(srcdir)/ChanContinuation.C 

chanContinuation_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
	$(top_builddir)/src-lapack/libnoxlapack.a \
	$(top_builddir)/test/utils/libnoxtestutils.a \
	$(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-lapack/liblocalapack.a \
	$(top_builddir)/test/lapack/LOCA_TestProblems/liblocatestproblems.a

chanContinuation_LDADD = $(TEST_LIBS)

chanTPContinuation_SOURCES = $(srcdir)/ChanTPContinuation.C 

chanTPContinuation_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanTPContinuation_LDADD = $(TEST_LIBS)

chanHomotopy_SOURCES = $(srcdir)/ChanHomotopy.C 

chanHomotopy_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanHomotopy_LDADD = $(TEST_LIBS)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)

