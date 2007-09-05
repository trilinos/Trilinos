# @HEADER
# ************************************************************************
#
#                           Intrepid Package
#                 Copyright (2007) Sandia Corporation
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
# Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
#                    Denis Ridzal (dridzal@sandia.gov).
#
# ************************************************************************
# @HEADER


# Tests

EXEEXT = .exe

noinst_PROGRAMS = test_01 test_02

test_01_SOURCES = test_01.cpp
test_01_LDADD = $(top_builddir)/src/libintrepid.a
test_01_CXXFLAGS = -I$(top_srcdir)/src/ -I$(top_srcdir)/src/local -I$(srcdir)

test_02_SOURCES = test_02.cpp
test_02_LDADD = $(top_builddir)/src/libintrepid.a
test_02_CXXFLAGS = -I$(top_srcdir)/src/ -I$(top_srcdir)/src/local -I$(srcdir)
