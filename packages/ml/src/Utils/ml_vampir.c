/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
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
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/

/* This file allows the use of the Vampir utility.
   
   Classes are general categories, e.g., Setup, Solve, Aggregation.   Within a
   class, you define functions and associated state variables (just integer
   handles).   (Functions don't necessarily have to correspond to C functions.)
   For instance, you might have an Aggregation class that contains functions
   "Phase 1" and "Phase 2/3 Cleanup."

*/

#include "ml_vampir.h"

#ifdef ML_VAMPIR
/* global variables for different classes */
int ml_vt_genhierarchy_class;

/* global variables for different states */
int /*ml_vt_reitzinger_state;*/
           ml_vt_aggregating_state,
           ml_vt_building_coarse_T_state,
           ml_vt_build_Pe_state,
           ml_vt_smooth_Pe_state,
           ml_vt_make_coarse_A_with_rap_state,
           ml_vt_reitzinger_cleanup_state;

void ML_Vampir_Setup()
{
   /* Define different classes for symbols */
   VT_classdef("Create Hierarchy",&ml_vt_genhierarchy_class);

   /* Define the instrumentation symbols */
  /*VT_funcdef("Reitzinger",ml_vt_genhierarchy_class,&ml_vt_reitzinger_state);*/
   VT_funcdef("Nodal Aggregation",ml_vt_genhierarchy_class,
                   &ml_vt_aggregating_state);
   VT_funcdef("Building coarse T",ml_vt_genhierarchy_class,
                   &ml_vt_building_coarse_T_state);
   VT_funcdef("Building Pe",ml_vt_genhierarchy_class,
           &ml_vt_build_Pe_state);
   VT_funcdef("Smoothing Pe",ml_vt_genhierarchy_class,
           &ml_vt_smooth_Pe_state);
   VT_funcdef("Building Coarse A",ml_vt_genhierarchy_class,
           &ml_vt_make_coarse_A_with_rap_state);
   VT_funcdef("Finishing up Reitzinger",ml_vt_genhierarchy_class,
           &ml_vt_reitzinger_cleanup_state);
}
#endif /*ifdef ML_VAMPIR*/
