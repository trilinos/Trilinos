/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/********************************************************************* */
/*          Decomposition with METIS                                  */
/********************************************************************* */

#ifndef __MLAGGPARMETIS__
#define __MLAGGPARMETIS__

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


  extern int ML_Aggregate_Set_ReqLocalCoarseSize( int, ML_Aggregate *ag,
	 				 int level,
					 int desired_aggre_per_proc );
  extern int ML_DecomposeGraph_BuildOffsets( int N_parts,
				    int offsets[],
				    int N_procs, USR_COMM );
  extern int ML_Aggregate_Set_OptimalNumberOfNodesPerAggregate( int optimal_value );
  extern int ML_Aggregate_CoarsenParMETIS( ML_Aggregate *ml_ag,
					   ML_Operator *Amatrix,
					   ML_Operator **Pmatrix,
					   ML_Comm *comm);
  extern int ML_Aggregate_Set_OptimalNumberOfNodesPerAggregate( int optimal_value );
  extern int ML_Aggregate_Get_OptimalNumberOfNodesPerAggregate( );

  extern int ML_CountNodesPerAggre(int Nrows, int GraphDecomposition[],
					int Naggre, int * NnodesPerAggre,
					USR_COMM Comm);

  extern int ML_BuildReorderedOffset( int starting_offset[],
				    int desired_aggre_per_proc, int Nprocs,
				    int nodes_per_aggre[], int Naggregates,
				    int reordered_offset[], int mypid );

  /* those are coded in ml_agg_METIS.c */
  extern int ML_Aggregate_Set_UseDropping(int i);

  extern int ML_Aggregate_Get_UseDropping();

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef __MLAGGPARMETIS__ */
