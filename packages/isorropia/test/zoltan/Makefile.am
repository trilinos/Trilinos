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

noinst_PROGRAMS = rebalance_zoltan

rebalance_zoltan_SOURCES = $(srcdir)/rebalance_zoltan.cpp

rebalance_zoltan_DEPENDENCIES = \
	$(top_builddir)/src/libisorropia.a \
	$(top_builddir)/utils/libisorropia_test_utils.a

rebalance_zoltan_CXXFLAGS = -I$(top_srcdir)/utils \
   $(ISORROPIA_INCLUDES)
rebalance_zoltan_LDADD = \
   -L$(top_builddir)/utils -lisorropia_test_utils $(ISORROPIA_LIBS)


