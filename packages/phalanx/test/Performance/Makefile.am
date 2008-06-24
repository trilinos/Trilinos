#@HEADER
# ************************************************************************
# 
#                     New_Package Example Package
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

include $(top_builddir)/Makefile.export.phalanx

# The following line helps the test harness recover from build errors.

all-local:

if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(PHALANX_TEST_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(PHALANX_TEST_INCLUDES))
else
TEST_LIBS     = $(PHALANX_TEST_LIBS)
TEST_INCLUDES = $(PHALANX_TEST_INCLUDES)
endif

AM_CPPFLAGS = -I$(srcdir)/../Utilities -I$(srcdir)/../Utilities/evaluators  $(TEST_INCLUDES)

EXEEXT = .exe

if BUILD_TESTS
noinst_PROGRAMS = performance_bracket_operator performance_allocators
endif

performance_bracket_operator_SOURCES = \
	$(srcdir)/Performance_BracketOperator.cpp 

performance_bracket_operator_DEPENDENCIES = \
	$(top_builddir)/src/libphalanx.a

performance_bracket_operator_LDADD = $(TEST_LIBS)

performance_allocators_SOURCES = \
	$(srcdir)/Performance_Allocators.cpp 

performance_allocators_DEPENDENCIES = \
	$(top_builddir)/test/Utilities/libphalanx_utilities.a \
	$(top_builddir)/src/libphalanx.a

performance_allocators_LDADD = $(TEST_LIBS)
