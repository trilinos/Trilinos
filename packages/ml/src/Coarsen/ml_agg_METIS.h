/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/********************************************************************* */
/*          Decomposition with METIS                                  */
/********************************************************************* */

#ifndef __MLAGGMETIS__
#define __MLAGGMETIS__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*MS*/
#define ML_AGGREGATE_OPTIONS_ID 13579

/* undefined will be a negative number */
#define ML_NUM_LOCAL_AGGREGATES     0
#define ML_NUM_GLOBAL_AGGREGATES    1
#define ML_NUM_NODES_PER_AGGREGATE  2

typedef struct ML_Aggregate_Options_Struct
{
  int id;
  int Naggregates_local;
  int Nnodes_per_aggregate;
  int Naggregates_global;
  int choice;
  int reordering_flag;
  int desired_aggre_per_proc; /* for ParMETIS */
} ML_Aggregate_Options;
/*ms*/

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

  extern int ML_Aggregate_Options_Defaults( ML_Aggregate_Options * pointer,
					    int NumLevels );

  extern int ML_Aggregate_Set_NodesPerAggr( ML *ml, ML_Aggregate *ag,
					    int level, int nodes_per_aggre );
  extern int ML_Aggregate_Set_LocalNumber( ML *ml, ML_Aggregate *ag,
					   int level, int Nlocal  );
  extern int ML_Aggregate_Set_GlobalNumber( ML *ml, ML_Aggregate *ag,
					    int level, int Nglobal  );
  extern int ML_Aggregate_Set_ReorderingFlag( ML *ml, ML_Aggregate *ag,
					      int level, int reordering_flag);
  extern int ML_Aggregate_CoarsenMETIS( ML_Aggregate *ml_ag,
					ML_Operator *Amatrix,
					ML_Operator **Pmatrix, ML_Comm *comm);
  extern int ML_DecomposeGraph_BuildOffsets( int N_parts,
					     int offsets[],
					     int N_procs,
					     USR_COMM comm);
  extern int ML_Aggregate_Set_OptimalNumberOfNodesPerAggregate( int optimal_value );
  extern int ML_Aggregate_Get_OptimalNumberOfNodesPerAggregate( );

  extern int ML_Aggregate_Set_UseDropping(int i);

  extern int ML_Aggregate_Get_UseDropping();


#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef __MLAGGMETIS__ */
