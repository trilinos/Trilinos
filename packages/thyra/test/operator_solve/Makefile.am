#@HEADER
# ************************************************************************
# 
#               Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms
#                 Copyright (2004) Sandia Corporation
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

EXEEXT = .exe

#
# Actual programs
#

TEST_LINEAR_OP_WITH_SOLVE = test_linear_op_with_solve

noinst_PROGRAMS = $(TEST_LINEAR_OP_WITH_SOLVE)

test_linear_op_with_solve_SOURCES = $(srcdir)/test_linear_op_with_solve.cpp
test_linear_op_with_solve_DEPENDENCIES = $(common_dependencies)
test_linear_op_with_solve_LDADD = $(common_ldadd)

#
# Common stuff
#

common_dependencies = \
  $(top_builddir)/src/libthyra.a

include $(top_builddir)/Makefile.export.thyra

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(THYRA_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(THYRA_LIBS))
else
EXPORT_INCLUDES = $(THYRA_INCLUDES)
EXPORT_LIBS = $(THYRA_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

# This is already added part of THYRA_LIBS and therefore automake does not need to add this again!
LIBS =
