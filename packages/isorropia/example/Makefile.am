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

UTIL_INCLUDE=-I$(top_srcdir)/utils
UTIL_LIB=$(top_builddir)/utils/libisorropia_test_utils.a

EXEEXT = .exe

SUBDIRS=geometric levelSchedule

if BUILD_EXAMPLES
PROGS = matrix_1 part_redist vert_weights graphedge_weights hgedge_weights
endif

noinst_PROGRAMS = $(PROGS)

matrix_1_SOURCES = $(srcdir)/matrix_1.cpp

matrix_1_DEPENDENCIES = \
	$(top_builddir)/src/libisorropia.a

matrix_1_CXXFLAGS = $(ISORROPIA_INCLUDES) $(UTIL_INCLUDE)
matrix_1_LDADD = $(ISORROPIA_LIBS) $(UTIL_LIB)

part_redist_SOURCES = $(srcdir)/part_redist.cpp

part_redist_DEPENDENCIES = \
        $(top_builddir)/src/libisorropia.a

part_redist_CXXFLAGS = $(ISORROPIA_INCLUDES) $(UTIL_INCLUDE)
part_redist_LDADD = $(ISORROPIA_LIBS) $(UTIL_LIB)

vert_weights_SOURCES = $(srcdir)/vert_weights.cpp

vert_weights_DEPENDENCIES = \
        $(top_builddir)/src/libisorropia.a

vert_weights_CXXFLAGS = $(ISORROPIA_INCLUDES) $(UTIL_INCLUDE)
vert_weights_LDADD = $(ISORROPIA_LIBS) $(UTIL_LIB)

graphedge_weights_SOURCES = $(srcdir)/graphedge_weights.cpp

graphedge_weights_DEPENDENCIES = \
        $(top_builddir)/src/libisorropia.a

graphedge_weights_CXXFLAGS = $(ISORROPIA_INCLUDES) $(UTIL_INCLUDE)
graphedge_weights_LDADD = $(ISORROPIA_LIBS) $(UTIL_LIB)

hgedge_weights_SOURCES = $(srcdir)/hgedge_weights.cpp

hgedge_weights_DEPENDENCIES = \
        $(top_builddir)/src/libisorropia.a

hgedge_weights_CXXFLAGS = $(ISORROPIA_INCLUDES) $(UTIL_INCLUDE)
hgedge_weights_LDADD = $(ISORROPIA_LIBS) $(UTIL_LIB)

