#@HEADER
# ************************************************************************
# 
#               Isorropia: Partitioning and Load Balancing Package
#                 Copyright (2006) Sandia Corporation
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
# ************************************************************************
#@HEADER

include $(top_builddir)/Makefile.export.isorropia

EXEEXT = .exe

if BUILD_EXAMPLES
PROGS = example_rcb
endif

UTIL_INCLUDE=-I$(top_srcdir)/utils
UTIL_LIB=$(top_builddir)/utils/libisorropia_test_utils.a

noinst_PROGRAMS = $(PROGS)

example_rcb_SOURCES = $(srcdir)/example_rcb.cpp 

example_rcb_DEPENDENCIES = \
	$(top_builddir)/utils/libisorropia_test_utils.a \
	$(top_builddir)/src/libisorropia.a

example_rcb_CXXFLAGS = $(ISORROPIA_INCLUDES) $(UTIL_INCLUDE)
example_rcb_LDADD = $(ISORROPIA_LIBS) $(UTIL_LIB)

all-local:
	@if test -e *.coords; then : ; else cp $(srcdir)/*.coords .; fi

