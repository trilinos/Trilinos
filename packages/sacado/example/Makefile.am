## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                           Sacado Package
#                 Copyright (2004) Sandia Corporation
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
# Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
# (etphipp@sandia.gov).
# 
# ************************************************************************
# @HEADER

EXEEXT = .exe

noinst_PROGRAMS = \
	ad_example \
	dfad_example \
	dfad_dfad_example \
	sfad_example \
	dmfad_example \
	trad_example \
	trad_dfad_example \
	dtaylor_example \
	RunValgrind

dfad_example_SOURCES = $(srcdir)/dfad_example.cpp
dfad_example_CXXFLAGS = \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
dfad_example_DEPENDENCIES = $(top_builddir)/src/libsacado.a
dfad_example_LDADD = -L$(top_builddir)/src -lsacado $(LIBS) $(FLIBS)

ad_example_SOURCES = $(srcdir)/ad_example.cpp
ad_example_CXXFLAGS = $(dfad_example_CXXFLAGS)
ad_example_DEPENDENCIES = $(dfad_example_DEPENDENCIES)
ad_example_LDADD = $(dfad_example_LDADD)

dfad_dfad_example_SOURCES = $(srcdir)/dfad_dfad_example.cpp
dfad_dfad_example_CXXFLAGS = $(dfad_example_CXXFLAGS)
dfad_dfad_example_DEPENDENCIES = $(dfad_example_DEPENDENCIES)
dfad_dfad_example_LDADD = $(dfad_example_LDADD)

sfad_example_SOURCES = $(srcdir)/sfad_example.cpp
sfad_example_CXXFLAGS =  $(dfad_example_CXXFLAGS)
sfad_example_DEPENDENCIES = $(dfad_example_DEPENDENCIES)
sfad_example_LDADD = $(dfad_example_LDADD)

dmfad_example_SOURCES = $(srcdir)/dmfad_example.cpp
dmfad_example_CXXFLAGS =  $(dfad_example_CXXFLAGS)
dmfad_example_DEPENDENCIES = $(dfad_example_DEPENDENCIES)
dmfad_example_LDADD = $(dfad_example_LDADD)

trad_example_SOURCES = $(srcdir)/trad_example.cpp
trad_example_CXXFLAGS = $(dfad_example_CXXFLAGS)
trad_example_DEPENDENCIES = $(dfad_example_DEPENDENCIES)
trad_example_LDADD = $(dfad_example_LDADD)

trad_dfad_example_SOURCES = $(srcdir)/trad_dfad_example.cpp
trad_dfad_example_CXXFLAGS = $(dfad_example_CXXFLAGS)
trad_dfad_example_DEPENDENCIES = $(dfad_example_DEPENDENCIES)
trad_dfad_example_LDADD = $(dfad_example_LDADD)

dtaylor_example_SOURCES = $(srcdir)/dtaylor_example.cpp
dtaylor_example_CXXFLAGS = $(dfad_example_CXXFLAGS)
dtaylor_example_DEPENDENCIES = $(dfad_example_DEPENDENCIES)
dtaylor_example_LDADD = $(dfad_example_LDADD)

RunValgrind_SOURCES = \
	$(srcdir)/RunValgrind \
	$(srcdir)/In_Suppressions

RunValgrind$(EXEEXT): $(srcdir)/RunValgrind
	cp  $(srcdir)/RunValgrind  RunValgrind$(EXEEXT)
	chmod +x RunValgrind$(EXEEXT)
	cp  $(srcdir)/In_Suppressions Suppressions
if BUILD_FEMEX
FEMDIR = 1DFEM
else
FEMDIR =
endif

SUBDIRS = $(FEMDIR)
