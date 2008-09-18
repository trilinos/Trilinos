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

if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(PHALANX_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(PHALANX_INCLUDES))
else
TEST_LIBS     = $(PHALANX_LIBS)
TEST_INCLUDES = $(PHALANX_INCLUDES)
endif

AM_CPPFLAGS = -I$(srcdir)/../Utilities $(TEST_INCLUDES)

EXEEXT = .exe

if BUILD_TESTS
noinst_PROGRAMS = data_container 
endif

data_container_SOURCES = \
	$(srcdir)/DataContainer.cpp 

data_container_DEPENDENCIES = \
	$(top_builddir)/src/libphalanx.a

data_container_LDADD = $(TEST_LIBS)
