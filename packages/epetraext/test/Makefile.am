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

if BUILD_TESTS

if BUILD_TRANSFORM

if BUILD_EXPERIMENTAL

TRANSFORM_SUBDIRS = Composite Permutation SingletonFilter View Overlap Transpose SymmRCM scripts

else

TRANSFORM_SUBDIRS = Composite Permutation SingletonFilter View Overlap Transpose scripts

endif

else

TRANSFORM_SUBDIRS =

endif

if BUILD_BTF

BTF_SUBDIRS = BTF

else

BTF_SUBDIRS =

endif

if BUILD_INOUT

INOUT_SUBDIRS = inout

else 

INOUT_SUBDIRS = 

endif 

if BUILD_ZOLTAN

ZOLTAN_SUBDIRS = Zoltan

else

ZOLTAN_SUBDIRS =

endif

if BUILD_COLORING

COLORING_SUBDIRS = MapColoring

else

COLORING_SUBDIRS = 

endif

SUBDIRS = MatrixMatrix $(TRANSFORM_SUBDIRS) $(BTF_SUBDIRS) $(INOUT_SUBDIRS) $(ZOLTAN_SUBDIRS) $(COLORING_SUBDIRS)

endif
