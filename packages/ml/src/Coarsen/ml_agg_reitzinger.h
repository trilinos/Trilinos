
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */      
  
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML_Operator structure                             */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __MLAGGREITZINGER__
#define __MLAGGREITZINGER__
#include "ml_defs.h"
#include "ml_mat_formats.h"
#include "ml_agg_genP.h"
#include "ml_op_utils.h"
/* ******************************************************************** */
/* ******************************************************************** */
/*      User Interface Proto-types                                      */
/* ******************************************************************** */
/* ******************************************************************** */

#ifdef __cplusplus
extern "C" {
#endif

extern int  ML_Gen_MGHierarchy_UsingReitzinger(ML *ml_edges, ML* ml_nodes, 
                                      int fine_level, int incr_or_decrease,
                                      ML_Aggregate *ag, ML_Operator *Tmat,
                                      ML_Operator *Tmat_trans,
                                      ML_Operator ***Tmat_array,
                                      ML_Operator ***Tmat_trans_array,
                                      int smooth_flag, double smooth_factor);
extern int ML_MGHierarchy_ReitzingerDestroy(int finest_level, 
                        int coarsest_level, ML_Operator ***Tmat_array,
                        ML_Operator ***Tmat_trans_array);

extern int ML_Gen_Hierarchy_ComplexMaxwell(ML *ml_edges,
					   ML_Operator **Tmat_array, 
					   ML_Operator **Tmat_trans_array,
					   ML **newml, ML_Operator *M);
#ifdef __cplusplus
}
#endif

#endif




