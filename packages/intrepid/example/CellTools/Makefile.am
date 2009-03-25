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

include $(top_builddir)/Makefile.export.intrepid

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(INTREPID_TEST_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(INTREPID_TEST_INCLUDES))
else
EXPORT_LIBS = $(INTREPID_TEST_LIBS)
EXPORT_INC_PATH = $(INTREPID_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

noinst_PROGRAMS = example_01 example_02 example_03 example_04

# Dependencies used when rebuilding examples without building everything from scratch.

example_01_DEPENDENCIES = \
	$(top_builddir)/src/libintrepid.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

example_01_SOURCES = example_01.cpp
example_01_LDADD = $(EXPORT_LIBS) 
example_01_CXXFLAGS = -I$(srcdir)

example_02_SOURCES = example_02.cpp
example_02_LDADD = $(EXPORT_LIBS) 
example_02_CXXFLAGS = -I$(srcdir)

example_03_SOURCES = example_03.cpp
example_03_LDADD = $(EXPORT_LIBS) 
example_03_CXXFLAGS = -I$(srcdir)

example_04_SOURCES = example_04.cpp
example_04_LDADD = $(EXPORT_LIBS) 
example_04_CXXFLAGS = -I$(srcdir)