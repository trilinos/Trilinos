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
#  CVS Information
#  $Source$
#  $Author$
#  $Date$
#  $Revision$
# ************************************************************************
#@HEADER

include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(NOX_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(NOX_INCLUDES))
else
TEST_LIBS     = $(NOX_LIBS)
TEST_INCLUDES = $(NOX_INCLUDES)
endif

AM_CPPFLAGS = $(TEST_INCLUDES)

EXEEXT = .exe

if BUILD_EPETRA_TESTS

if BUILD_NOX_EPETRAEXT_SUPPORT
noinst_PROGRAMS = FiniteDifferenceColoring
endif

endif

FiniteDifferenceColoring_SOURCES = \
        $(srcdir)/Example.C \
	$(srcdir)/FiniteElementProblem.C \
        $(srcdir)/FiniteElementProblem.H \
	$(srcdir)/Basis.C \
        $(srcdir)/Basis.H \
	$(srcdir)/Problem_Interface.C \
        $(srcdir)/Problem_Interface.H 

FiniteDifferenceColoring_DEPENDENCIES = $(top_builddir)/src/libnox.a \
        $(top_builddir)/src-epetra/libnoxepetra.a

FiniteDifferenceColoring_LDADD = $(TEST_LIBS)
