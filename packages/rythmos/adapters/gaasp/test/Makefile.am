#@HEADER
# ***********************************************************************
#
#                           Rythmos Package
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
# Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
#
# ***********************************************************************
#@HEADER

EXEEXT = .exe

include $(top_builddir)/Makefile.export.rythmos

bin_PROGRAMS = test_gaasp_adapters

#
# Actual programs
#

#test_gaasp_adapters_INCLUDES =
test_gaasp_adapters_SOURCES =\
	$(srcdir)/test_gaasp_adapters.cpp\
	$(srcdir)/Rythmos_GAASPErrorEstimator_UnitTest.cpp\
	$(srcdir)/../example/LorenzModel.cpp\
test_gaasp_adapters_DEPENDENCIES = $(top_builddir)/src/librythmos.a 
test_gaasp_adapters_LDADD = $(RYTHMOS_GAASP_LIBS) $(GAASP_LIBS) $(RYTHMOS_LIBS) 

AM_CPPFLAGS = $(RYTHMOS_GAASP_INCLUDES) $(GAASP_INCLUDES) $(RYTHMOS_INCLUDES) 

# This is already added part of RYTHMOS_LIBS and therefore automake does not need to add this again!
LIBS =
