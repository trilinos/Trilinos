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

noinst_PROGRAMS = copyfiles utest

utest_SOURCES = $(srcdir)/utest.cpp $(srcdir)/utst_serial.cpp $(srcdir)/utst_serial.hpp

copyfiles_SOURCES = $(top_srcdir)/test/utest/src_utest_dist1

utest_DEPENDENCIES = \
	$(top_builddir)/src/libisorropia.a \
	$(top_builddir)/utils/libisorropia_test_utils.a

utest_CXXFLAGS = -I$(top_srcdir)/utils $(ISORROPIA_INCLUDES) 
utest_LDADD = -L$(top_builddir)/utils -lisorropia_test_utils \
     $(ISORROPIA_LIBS) 

EXTRA_DIST = \
	src_utest_dist1 \
        src_utest_dist2 \
        src_utest_dist3

copyfiles$(EXEEXT):$(copyfiles_SOURCES)
	cp $(top_srcdir)/test/utest/src_utest_dist1 $(top_builddir)/test/utest/utest_dist1
	cp $(top_srcdir)/test/utest/src_utest_dist2 $(top_builddir)/test/utest/utest_dist2
	cp $(top_srcdir)/test/utest/src_utest_dist3 $(top_builddir)/test/utest/utest_dist3

