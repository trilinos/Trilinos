#@HEADER
# ************************************************************************
# 
#        Phalanx: A Partial Differential Equation Field Evaluation 
#       Kernel for Flexible Management of Complex Dependency Chains
#                  Copyright (2008) Sandia Corporation
# 
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
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
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
# National Laboratories.
# 
# ************************************************************************
#@HEADER

include $(top_builddir)/Makefile.export.phalanx

# The following line helps the test harness recover from build errors.

all-local:

## Do not use PHALANX_TEST_<>, but use PHALANX_<>.  Don't want to pull in
## the test/Utilities.
if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(PHALANX_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(PHALANX_INCLUDES))
else
TEST_LIBS     = $(PHALANX_LIBS)
TEST_INCLUDES = $(PHALANX_INCLUDES)
endif

## Do not point to utilities directory.  Use special algebraic type file
## in this directory.
AM_CPPFLAGS = $(TEST_INCLUDES)

EXEEXT = .exe

if BUILD_TESTS
noinst_PROGRAMS = performance_algebraic_types performance_algebraic_types_multi
endif

performance_algebraic_types_SOURCES = \
	$(srcdir)/Performance_AlgebraicTypes.cpp \
	$(srcdir)/Special_AlgebraicTypes.hpp

performance_algebraic_types_DEPENDENCIES = \
	$(top_builddir)/src/libphalanx.a

performance_algebraic_types_LDADD = $(TEST_LIBS)

performance_algebraic_types_multi_SOURCES = \
	$(srcdir)/Performance_AlgebraicTypes_Multi.cpp \
	$(srcdir)/Special_AlgebraicTypes.hpp

performance_algebraic_types_multi_DEPENDENCIES = \
	$(top_builddir)/src/libphalanx.a

performance_algebraic_types_multi_LDADD = $(TEST_LIBS)
