
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
#include "ml_common.h"
#include "ml_defs.h"
#include "ml_mat_formats.h"
#include "ml_agg_genP.h"
#include "ml_op_utils.h"
#include "ml_operator_blockmat.h"
#include "ml_utils.h"
/* ******************************************************************** */
/* ******************************************************************** */
/*      User Interface Proto-types                                      */
/* ******************************************************************** */
/* ******************************************************************** */
struct ml_linked_list {
  struct ml_linked_list *next;
  int duplicate_row;
};

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

extern int  ML_Gen_MGHierarchy_UsingReitzinger(ML *ml_edges, ML* ml_nodes, 
                                      int fine_level, int incr_or_decrease,
                                      ML_Aggregate *ag, ML_Operator *Tmat,
                                      ML_Operator *Tmat_trans,
                                      ML_Operator ***Tmat_array,
                                      ML_Operator ***Tmat_trans_array,
                                      int smooth_flag, double smooth_factor);
extern int ML_MGHierarchy_ReitzingerDestroy(int finest_level, 
                        ML_Operator ***Tmat_array,
                        ML_Operator ***Tmat_trans_array);

extern int ML_Gen_Hierarchy_ComplexMaxwell(ML *ml_edges,
					   ML **newml, ML_Operator *M);

extern int ML_Gen_SmoothPnodal(ML *ml,int level, int clevel, void *data,
			       double smoothP_damping_factor,
			       ML_Operator *SPn_mat);

 
extern int ml_comp_Pe_entries(int coef_cols[], double coef_values[], 
			      int coef_count, int leftagg, 
			      struct ml_linked_list **Trecorder,
			      int *Trowcount, int *Tnzcount,
			      struct ML_CSR_MSRdata *Tcoarse,
			      int *Pnzcount, int Pe_columns[], 
			      double Pe_values[]);

extern int ml_record_entry(struct ml_linked_list **Trecorder,int lower, 
			int therow);
extern int ml_dup_entry(int node1, int node2, struct ml_linked_list **Trecorder,
	int Tcols[], int Trowptr[], int *lower, int *upper, 
	int *duplicate_row);

extern int ml_clean_Trecorder(struct ml_linked_list ***Trecorder ,int N);

extern int ml_leastsq_edge_interp(ML_Operator *Pn_mat, ML_Operator *SPn_mat, 
			   ML_Operator *Tfine_mat, ML_Operator *Tcoarse_mat, 
			   ML_Operator *Pe_mat, int);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif




