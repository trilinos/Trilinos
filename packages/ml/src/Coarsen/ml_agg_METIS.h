/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */   
/* ******************************************************************** */

/********************************************************************* */
/*          Decomposition with METIS                                  */
/********************************************************************* */

#ifndef __MLAGGMETIS__
#define __MLAGGMETIS__

/*MS*/
#define ML_AGGREGATE_OPTIONS_ID 13579

typedef struct ML_Aggregate_Options_Struct
{
  int id;
  int Naggregates;
  int Nnodes_per_aggregate;
  int local_or_global;
  int reordering_flag;
} ML_Aggregate_Options;
/*ms*/

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif
  
  extern int ML_Aggregate_Set_NodesPerAggr( ML *ml, ML_Aggregate *ag, 
					    int level, int nodes_per_aggre );
  extern int ML_Aggregate_Set_LocalNumber( ML *ml, ML_Aggregate *ag, 
					   int level, int Nlocal  );
  extern int ML_Aggregate_Set_ReorderingFlag( ML *ml, ML_Aggregate *ag, 
					      int level, int reordering_flag);
  extern int ML_Aggregate_CoarsenMETIS( ML_Aggregate *ml_ag,
					ML_Operator *Amatrix, 
					ML_Operator **Pmatrix, ML_Comm *comm);
  extern int ML_DecomposeGraph_BuildOffsets( int N_parts,
					     int offsets[],
					     int N_procs );
  extern int ML_Aggregates_CheckAggregates( int Naggregates, int N_rows,
					    int graph_decomposition[],
					    int mypid);
  

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef __MLAGGMETIS__ */
