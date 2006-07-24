# @HEADER
# ************************************************************************
# 
#                MLMEX: MATLAB(Mex) Interface for ML/MLAPI
#                   Copyright (2006) Sandia Corporation
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

# Include export Makefiles
include  $(top_builddir)/Makefile.export

# I need extra_libs because of MLAPI::Gallery, which calls Galeri.
ML_MLMEX_INCLUDES = $(ML_INCLUDES) $(ML_EXTRA_INCLUDES) $(ML_MATLAB_INCLUDES)
ML_MLMEX_LIBS     = $(ML_LIBS) $(ML_EXTRA_LIBS)

# The MEXEXT_VALUE is set in AZ_MATLAB_ROOT from the mexext script, the only portable
# way of figuring out what my mex system extension is.  Of course, this means cross-compiling
# probably won't work.  C'est la vie.
EXEEXT = .$(MEXEXT_VALUE)

if MATLAB_USE
bin_PROGRAMS = mlmex
bin_SCRIPTS = ml.m

mlmex_SOURCES      = mlmex.cpp
mlmex_DEPENDENCIES = $(top_builddir)/src/libml.a
mlmex_CXXFLAGS     = ${ML_MLMEX_INCLUDES} 
mlmex_LDADD        = ${ML_MLMEX_LIBS} 
mlmex_LINK         = $(MEX) -cxx -o $@
endif
