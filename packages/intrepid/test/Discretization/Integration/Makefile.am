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

noinst_PROGRAMS = \
test_01 \
test_02 \
test_03 \
test_04 \
test_05 \
test_06 \
test_07 \
test_08 \
test_09 \
test_10 \
Dummy

# Dependencies used when rebuilding tests without building everything from scratch.
test_01_DEPENDENCIES = \
	$(top_builddir)/src/libintrepid.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../shards/src/libshards.a
test_02_DEPENDENCIES = \
	$(top_builddir)/src/libintrepid.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../shards/src/libshards.a
test_03_DEPENDENCIES = \
	$(top_builddir)/src/libintrepid.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../shards/src/libshards.a
test_04_DEPENDENCIES = \
	$(top_builddir)/src/libintrepid.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../shards/src/libshards.a
test_05_DEPENDENCIES = \
	$(top_builddir)/src/libintrepid.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../shards/src/libshards.a
test_06_DEPENDENCIES = \
	$(top_builddir)/src/libintrepid.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../shards/src/libshards.a
test_07_DEPENDENCIES = \
	$(top_builddir)/src/libintrepid.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../shards/src/libshards.a
test_08_DEPENDENCIES = \
	$(top_builddir)/src/libintrepid.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../shards/src/libshards.a
test_09_DEPENDENCIES = \
	$(top_builddir)/src/libintrepid.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../shards/src/libshards.a
test_10_DEPENDENCIES = \
	$(top_builddir)/src/libintrepid.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../shards/src/libshards.a

test_01_SOURCES = test_01.cpp
test_01_LDADD = $(EXPORT_LIBS)
test_01_CXXFLAGS = -I$(srcdir)

test_02_SOURCES = test_02.cpp
test_02_LDADD = $(EXPORT_LIBS)
test_02_CXXFLAGS = -I$(srcdir)

test_03_SOURCES = test_03.cpp
test_03_LDADD = $(EXPORT_LIBS)
test_03_CXXFLAGS = -I$(srcdir)

test_04_SOURCES = test_04.cpp
test_04_LDADD = $(EXPORT_LIBS)
test_04_CXXFLAGS = -I$(srcdir)

test_05_SOURCES = test_05.cpp
test_05_LDADD = $(EXPORT_LIBS)
test_05_CXXFLAGS = -I$(srcdir)

test_06_SOURCES = test_06.cpp
test_06_LDADD = $(EXPORT_LIBS)
test_06_CXXFLAGS = -I$(srcdir)

test_07_SOURCES = test_07.cpp
test_07_LDADD = $(EXPORT_LIBS)
test_07_CXXFLAGS = -I$(srcdir)

test_08_SOURCES = test_08.cpp
test_08_LDADD = $(EXPORT_LIBS)
test_08_CXXFLAGS = -I$(srcdir)

test_09_SOURCES = test_09.cpp
test_09_LDADD = $(EXPORT_LIBS)
test_09_CXXFLAGS = -I$(srcdir)

test_10_SOURCES = test_10.cpp
test_10_LDADD = $(EXPORT_LIBS)
test_10_CXXFLAGS = -I$(srcdir)

Dummy_SOURCES = Dummy
Dummy$(EXEEXT): $(srcdir)/Dummy
	if [ $(srcdir) != '.' ]; then mkdir data; cp -r $(srcdir)/./data/*.dat data; fi

EXTRA_DIST = \
data/HEX_integrals.dat       data/TET_integrals.dat \
data/TRIPRISM_integrals.dat  data/EDGE_integrals.dat \
data/QUAD_integrals.dat      data/TRI_integrals.dat
