## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                           Sacado Package
#                 Copyright (2006) Sandia Corporation
# 
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
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

all-local:
	@echo ""
	@echo "Trilinos package sacado subdirectory example built successfully."
	@echo ""

noinst_PROGRAMS = \
	ad_example \
	dfad_example \
	dfad_dfad_example \
	sfad_example \
	dmfad_example \
	trad_example \
	trad_dfad_example \
	taylor_example \
	dfad_taylor_example

dfad_example_SOURCES = $(srcdir)/dfad_example.cpp
dfad_example_CXXFLAGS = \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/src/mpl \
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

taylor_example_SOURCES = $(srcdir)/taylor_example.cpp
taylor_example_CXXFLAGS = $(dfad_example_CXXFLAGS)
taylor_example_DEPENDENCIES = $(dfad_example_DEPENDENCIES)
taylor_example_LDADD = $(dfad_example_LDADD)

dfad_taylor_example_SOURCES = $(srcdir)/dfad_taylor_example.cpp
dfad_taylor_example_CXXFLAGS = $(dfad_example_CXXFLAGS)
dfad_taylor_example_DEPENDENCIES = $(dfad_example_DEPENDENCIES)
dfad_taylor_example_LDADD = $(dfad_example_LDADD)

include FEApp/MakefileFEAppStub.am
EXTRA_DIST = $(FEAPP_H) $(FEAPP)
