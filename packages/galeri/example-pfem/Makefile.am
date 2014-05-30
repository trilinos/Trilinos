
# The following line helps the test harness recover from build errors.DER
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
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @HEADER

all-local:
	@echo ""
	@echo "Trilinos package galeri subdirectory example-pfem built successfully."
	@echo ""

include $(top_builddir)/Makefile.export.galeri

EXEEXT = .exe

noinst_PROGRAMS = ex1 ex2 ex3

ex1_SOURCES      = $(srcdir)/ex1.cpp
ex1_DEPENDENCIES = $(top_builddir)/src-pfem/libpfem.a
ex1_CXXFLAGS     = $(GALERI_INCLUDES)
ex1_LDADD        = $(GALERI_LIBS) 

ex2_SOURCES      = $(srcdir)/ex2.cpp
ex2_DEPENDENCIES = $(top_builddir)/src-pfem/libpfem.a
ex2_CXXFLAGS     = $(GALERI_INCLUDES)
ex2_LDADD        = $(GALERI_LIBS) 

ex3_SOURCES      = $(srcdir)/ex3.cpp
ex3_DEPENDENCIES = $(top_builddir)/src-pfem/libpfem.a
ex3_CXXFLAGS     = $(GALERI_INCLUDES)
ex3_LDADD        = $(GALERI_LIBS) 
