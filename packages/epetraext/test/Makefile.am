#@HEADER
# ************************************************************************
# 
#     EpetraExt: Epetra Extended - Linear Algebra Services Package 
#                 Copyright (2001) Sandia Corporation
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

# The following line helps the test harness recover from build errors.

all-local:

if BUILD_TESTS

if BUILD_TRANSFORM

if BUILD_GRAPH_REORDERINGS
if BUILD_EXPERIMENTAL
COMPOSITE_SUBDIR = Composite
SYMMRCM_SUBDIR = SymmRCM
endif
endif

PERMUTATION_SUBDIR = Permutation
SINGLETONFILTER_SUBDIR = SingletonFilter
VIEW_SUBDIR = View
OVERLAP_SUBDIR = Overlap
if BUILD_EPETRAEXT_TRIUTILS
TRANSPOSE_SUBDIR = Transpose
endif

endif

TRANSFORM_SUBDIRS = $(COMPOSITE_SUBDIR) $(PERMUTATION_SUBDIR) \
	$(SINGLETONFILTER_SUBDIR) $(VIEW_SUBDIR) $(OVERLAP_SUBDIR) \
	$(TRANSPOSE_SUBDIR) $(SYMMRCM_SUBDIR) 

if BUILD_BTF

BTF_SUBDIRS = BTF

endif

if BUILD_INOUT
if BUILD_EPETRAEXT_TRIUTILS
INOUT_SUBDIRS = inout
endif
endif 

if BUILD_ZOLTAN

ZOLTAN_SUBDIRS = Zoltan

endif

if BUILD_COLORING

COLORING_SUBDIRS = MapColoring

endif

if BUILD_MATLAB
if BUILD_EPETRAEXT_TRIUTILS
MATLAB_SUBDIRS = Matlab
endif
endif

SUBDIRS = \
  MatrixMatrix \
  $(TRANSFORM_SUBDIRS) $(BTF_SUBDIRS) \
  $(INOUT_SUBDIRS) $(ZOLTAN_SUBDIRS) $(COLORING_SUBDIRS) \
  $(MATLAB_SUBDIRS)

endif 

EXTRA_DIST = definition

