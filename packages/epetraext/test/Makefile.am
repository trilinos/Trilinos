#@HEADER
# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                Copyright (2001) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#   
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#   
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
# 
# ************************************************************************
#@HEADER

if BUILD_TESTS

if BUILD_TRANSFORM

if BUILD_EXPERIMENTAL

TRANSFORM_SUBDIRS = Composite Permutation SingletonFilter View Overlap SymmRCM scripts

else

TRANSFORM_SUBDIRS = Composite Permutation SingletonFilter View Overlap scripts

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
