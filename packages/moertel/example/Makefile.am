#@HEADER
# ************************************************************************
# 
#                          Moertel FE Package 
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
# Questions? Contact Michael Gee (mwgee@sandia.gov) 
# 
# ************************************************************************
#@HEADER

SUBDIRS =

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package moertel subdirectory example built successfully."
	@echo ""

include $(top_builddir)/Makefile.export.moertel
include $(top_builddir)/../galeri/Makefile.export.galeri

EXEEXT = .exe

noinst_PROGRAMS = TwoSquares Hex_3D test1

EXTRA_DIST = TwoSquares.grid Hex_3D.grid TwoSquares.m

TwoSquares_SOURCES      = $(srcdir)/TwoSquares.cpp
TwoSquares_DEPENDENCIES = $(top_builddir)/src/libmoertel.a
TwoSquares_CXXFLAGS     = $(MOERTEL_INCLUDES) $(GALERI_INCLUDES)
TwoSquares_LDADD        = $(MOERTEL_LIBS) $(GALERI_LIBS)

Hex_3D_SOURCES      = $(srcdir)/Hex_3D.cpp
Hex_3D_DEPENDENCIES = $(top_builddir)/src/libmoertel.a
Hex_3D_CXXFLAGS     = $(MOERTEL_INCLUDES) $(GALERI_INCLUDES)
Hex_3D_LDADD        = $(MOERTEL_LIBS) $(GALERI_LIBS)

test1_SOURCES      = $(srcdir)/test1.cpp
test1_DEPENDENCIES = $(top_builddir)/src/libmoertel.a
test1_CXXFLAGS     = $(MOERTEL_INCLUDES)
test1_LDADD        = $(MOERTEL_LIBS)
