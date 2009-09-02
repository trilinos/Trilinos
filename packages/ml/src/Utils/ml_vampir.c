/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
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
