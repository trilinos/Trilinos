# @HEADER
# ************************************************************************
#
#                  Galeri Matrix Generation Package
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package galeri subdirectory example-fem built successfully."
	@echo ""

include $(top_builddir)/Makefile.export.galeri

EXEEXT = .exe

noinst_PROGRAMS = Grid Laplacian3D AdvDiff2D

Grid_SOURCES      = $(srcdir)/Grid.cpp
Grid_DEPENDENCIES = $(top_builddir)/src/libgaleri.a
Grid_CXXFLAGS     = $(GALERI_INCLUDES)
Grid_LDADD        = $(GALERI_LIBS)

Laplacian3D_SOURCES      = $(srcdir)/Laplacian3D.cpp
Laplacian3D_DEPENDENCIES = $(top_builddir)/src/libgaleri.a
Laplacian3D_CXXFLAGS     = $(GALERI_INCLUDES) 
Laplacian3D_LDADD        = $(GALERI_LIBS)

AdvDiff2D_SOURCES      = $(srcdir)/AdvDiff2D.cpp
AdvDiff2D_DEPENDENCIES = $(top_builddir)/src/libgaleri.a
AdvDiff2D_CXXFLAGS     = $(GALERI_INCLUDES)
AdvDiff2D_LDADD        = $(GALERI_LIBS)
