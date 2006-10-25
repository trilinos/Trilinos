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

# Include NOX export Makefile for correct paths
include $(top_builddir)/../teuchos/Makefile.export.teuchos

EXEEXT = .exe

noinst_PROGRAMS = \
	dfad_nonlinear \
	dfad_mult

dfad_nonlinear_SOURCES = dfad_nonlinear.cpp 
dfad_nonlinear_CXXFLAGS = -O2 \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src \
	-I$(top_srcdir)/test/utils \
	-I$(top_srcdir)/Fad \
        $(TEUCHOS_INCLUDES)
dfad_nonlinear_DEPENDENCIES = \
        $(top_builddir)/src/libsacado.a \
	$(top_builddir)/test/utils/libsacadoutils.a
dfad_nonlinear_LDADD = \
	-L$(top_builddir)/src -lsacado \
	-L$(top_builddir)/test/utils -lsacadoutils \
	$(TEUCHOS_LIBS) \
	$(LIBS) $(FLIBS)

dfad_mult_SOURCES = dfad_mult.cpp
dfad_mult_CXXFLAGS = $(dfad_nonlinear_CXXFLAGS)
dfad_mult_DEPENDENCIES = $(dfad_nonlinear_DEPENDENCIES)
dfad_mult_LDADD = $(dfad_nonlinear_LDADD)
