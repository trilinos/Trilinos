#@HEADER
# ************************************************************************
# 
#                 ML
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

EXEEXT = .exe

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
ml_common_include_locations = $(ML_EXPORT_INCLUDES)
ml_libs =       $(ML_EXPORT_LIBS)     \
                $(LIBS)               \
                $(LAPACK_LIBS)        \
                $(BLAS_LIBS)          \
                $(FLIBS)

noinst_PROGRAMS = VerySimple

VerySimple_SOURCES = $(srcdir)/VerySimple.c

VerySimple_DEPENDENCIES = \
	$(top_builddir)/src/libml.a 

mllink = $(CXX) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

VerySimple_CPPFLAGS = ${ml_common_include_locations}
VerySimple_CXXFLAGS = ${ml_common_include_locations}
VerySimple_LDADD = ${ml_libs}
VerySimple_LINK = ${mllink}
