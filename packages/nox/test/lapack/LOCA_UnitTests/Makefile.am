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

# Include NOX export makefile
include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(LOCA_LAPACK_TEST_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(LOCA_LAPACK_TEST_INCLUDES))
else
TEST_LIBS     = $(LOCA_LAPACK_TEST_LIBS)
TEST_INCLUDES = $(LOCA_LAPACK_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(TEST_INCLUDES)

if BUILD_LOCAANASAZI

ANASAZIJACOBIANINVERSE = anasaziJacobianInverse_exe

else

ANASAZIJACOBIANINVERSE = 

endif

BORDEREDSOLVENC = borderedSolveNonContiguous_exe
BORDEREDSOLVEC = borderedSolveContiguous_exe

noinst_PROGRAMS = \
	$(ANASAZIJACOBIANINVERSE) \
	$(BORDEREDSOLVENC) \
	$(BORDEREDSOLVEC)

anasaziJacobianInverse_exe_SOURCES = \
	$(srcdir)/AnasaziJacobianInverse.C 

anasaziJacobianInverse_exe_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
	$(top_builddir)/src-lapack/libnoxlapack.a \
	$(top_builddir)/test/utils/libnoxtestutils.a \
	$(top_builddir)/src-loca/src/libloca.la \
        $(top_builddir)/src-loca/src-lapack/liblocalapack.la \
	$(top_builddir)/test/lapack/LOCA_TestProblems/liblocatestproblems.la

anasaziJacobianInverse_exe_LDADD = $(TEST_LIBS)

borderedSolveNonContiguous_exe_SOURCES = \
	$(srcdir)/BorderedSolveNonContiguous.C

borderedSolveNonContiguous_exe_DEPENDENCIES = \
	$(anasaziJacobianInverse_exe_DEPENDENCIES)

borderedSolveNonContiguous_exe_LDADD = $(TEST_LIBS)

borderedSolveContiguous_exe_SOURCES = \
	$(srcdir)/BorderedSolveContiguous.C

borderedSolveContiguous_exe_DEPENDENCIES = \
	$(anasaziJacobianInverse_exe_DEPENDENCIES)

borderedSolveContiguous_exe_LDADD = $(TEST_LIBS)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)

# ------------------------------------------------------------------------
# Include global libtool definitions
# ------------------------------------------------------------------------
include $(top_srcdir)/config/libtool_macros.am
