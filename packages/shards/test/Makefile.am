#@HEADER
# ************************************************************************
# 
#                 Shards:  Shared Discretization Tools
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

include $(top_builddir)/Makefile.export.shards

# The following line helps the test harness recover from build errors.

all-local:

if USING_GNUMAKE
TEST_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(SHARDS_LIBS))
TEST_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(SHARDS_INCLUDES))
else
TEST_LIBS     = $(SHARDS_LIBS)
TEST_INCLUDES = $(SHARDS_INCLUDES)
endif

AM_CPPFLAGS = $(TEST_INCLUDES)

EXEEXT = .exe

## if BUILD_TESTS
bin_PROGRAMS = test_all
## endif

test_all_SOURCES = \
	$(srcdir)/Test_Main.cpp	\
	$(srcdir)/Test_Shards_Array.cpp 

test_all_DEPENDENCIES = \
	$(top_builddir)/src/libshards.a

test_all_LDADD = $(TEST_LIBS)
