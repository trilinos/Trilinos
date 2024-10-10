/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/********************************************************************* */
/*          Decomposition with Zoltan                                 */
/********************************************************************* */

#ifndef ML_AGG_ZOLTAN_H
#define ML_AGG_ZOLTAN_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_comm.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

  extern int ML_Aggregate_CoarsenParZoltan(ML_Aggregate *ml_ag,
					   ML_Operator *Amatrix,
					   ML_Operator **Pmatrix,
					   ML_Comm *comm);
  extern int ML_Aggregate_CoarsenZoltan(ML_Aggregate *ml_ag,
                   ML_Operator *Amatrix, ML_Operator **Pmatrix, ML_Comm *comm);
  /* those are coded in ml_agg_METIS.c */
  extern int ML_Aggregate_Set_UseDropping(int i);

  extern int ML_Aggregate_Get_UseDropping();

  extern int ML_DecomposeGraph_with_Zoltan(ML_Operator *Amatrix,
				  int N_parts,
				  int graph_decomposition[],
				  double bdry_nodes[],
				  double [], double [], double [],
				  int,int,int,int,int,int);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef ML_AGG_ZOLTAN_H */
