/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* data structure to hold aggregation information                            */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL)                                       */
/* Date          : August, 1999                                              */
/* ************************************************************************* */

#ifndef __MLAGGH__
#define __MLAGGH__

#include <stdio.h>
#include <stdlib.h>

/* ************************************************************************* */
/* definition of the aggregate structure                                     */
/* ------------------------------------------------------------------------- */

typedef struct ML_Aggregate_Struct
{
   int    ML_id;
   double print_flag;
   int    max_coarse_size;             /* maximum size of coarsest grid */
   int    ordering;                    /* natural, random, graph        */
   int    min_nodes_per_aggregate;     /* aggregate size control        */
   int    max_neigh_already_selected;  /* complexity control            */
   double threshold;                   /* for pruning matrix            */ 
   double curr_threshold;              /* adjusted for levels           */
   double drop_tol_for_smoothing;      /* self-explanatory              */
   int    attach_scheme;               /* aggregate shape control       */
   int    spectral_radius_scheme;      /* way to compute approx maxeigen*/
   double smoothP_damping_factor;      /* for prolongator smoother      */
   int    smoothP_type;                /* point, block                  */
   int    coarsen_scheme;              /* Uncoupled, Coupled, MIS       */
   int    num_PDE_eqns;                /* block size                    */
   int    nullspace_dim;               /* self-explanatory              */
   double *nullspace_vect;             /* for null space vectors        */
   int    nullspace_corrupted;         /* indicates whether fine grid   */
                                       /* nullspace has been overwritten*/
   int    *aggr_count;                 /* no. aggregates at each level  */
   int    **aggr_info;                 /* node to aggregate map         */
   int    max_levels;                  /* maximum number of levels      */
   int    begin_level;                 /* finest grid level             */
   int    cur_level;                   /* temporary variable            */
   double operator_complexity;         /* sum of nnz for all A's        */
   double fine_complexity;             /* nnz of the finest A           */
   int    nvblocks;                    /* for variable blocks (finest)  */
   int    *vblock_info;                /* for variable blocks (finest)  */
   int    keep_P_tentative;            /* keeping tentative prolongator */
   void   *P_tentative;                /* so it can be reused later.    */
   int    smooth_existing_P_tentative; /* already have P tent, don't create it*/
} ML_Aggregate;

/* ************************************************************************* */
/* other ML include files                                               */
/* ------------------------------------------------------------------------- */

#include "ml_defs.h"
#include "ml_comm.h"
#include "ml_memory.h"
#include "ml_operator.h"

/* ************************************************************************* */
/* linked list structures for holding free node information             */
/* ------------------------------------------------------------------------- */

typedef struct ML_Node_Struct
{
   int    node_id;
   struct ML_Node_Struct *next;
} ML_Node;

/* ************************************************************************* */
/* definition of the structure for holding aggregate information        */
/* ------------------------------------------------------------------------- */

typedef struct ML_SuperNode_Struct
{
   int    length;
   int    maxlength;
   int    index;
   int    *list;
   struct ML_SuperNode_Struct *next;

} ML_SuperNode;

/* ************************************************************************* */
/* definition of the structure for holding communication information in */
/* aggregation procedure                                                */
/* ------------------------------------------------------------------------- */

typedef struct ML_Aggregate_Comm_Struct
{
   int     N_send_neighbors, N_recv_neighbors;
   int     local_nrows;
   int     *send_neighbors, *recv_neighbors;
   int     *send_leng, *recv_leng;
   int     *send_list;
   ML_Comm *comm;

} ML_Aggregate_Comm;

/* ************************************************************************* */
/* ************************************************************************* */
/* functions to manipulate the aggregate data structure                      */
/* ------------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" 
{
#endif

/* ------------------------------------------------------------------------- */
/* constructor/destructor and level control                                  */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Create( ML_Aggregate ** );
int ML_Aggregate_Destroy( ML_Aggregate ** );

int ML_Aggregate_Set_OutputLevel( ML_Aggregate *, double level );
int ML_Aggregate_Set_Reuse(ML_Aggregate *ag);

int ML_Aggregate_Set_MaxLevels( ML_Aggregate *, int level );
int ML_Aggregate_Set_CurrentLevel( ML_Aggregate *, int level );
int ML_Aggregate_Set_StartLevel( ML_Aggregate *, int level );

/* ------------------------------------------------------------------------- */
/* aggregate size and shape control                                          */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_MinNodesPerAggregate(ML_Aggregate *,int n);
int ML_Aggregate_Set_AttachScheme_MaxLink( ML_Aggregate * );
int ML_Aggregate_Set_AttachScheme_MinRank( ML_Aggregate * );

/* ------------------------------------------------------------------------- */
/* aggregation node traversal control                                        */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_NaturalOrdering( ML_Aggregate * );
int ML_Aggregate_Set_RandomOrdering( ML_Aggregate * );
int ML_Aggregate_Set_GraphOrdering( ML_Aggregate * );

/* ------------------------------------------------------------------------- */
/* control when to stop coarsening                                           */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_MaxCoarseSize( ML_Aggregate *, int size );

/* ------------------------------------------------------------------------- */
/* different parallel coarsening schemes                                     */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenScheme_Uncoupled( ML_Aggregate * );
int ML_Aggregate_Set_CoarsenScheme_Coupled( ML_Aggregate * );
int ML_Aggregate_Set_CoarsenScheme_MIS( ML_Aggregate *ag  );
int ML_Aggregate_Set_CoarsenScheme_DD( ML_Aggregate *ag  );
int ML_Aggregate_Set_CoarsenScheme_UncoupledMIS( ML_Aggregate *ag  );
int ML_Aggregate_Set_CoarsenScheme_UncoupledCoupled( ML_Aggregate *ag  );
int ML_Aggregate_Phase2_3_Cleanup(ML_Aggregate *ml_ag, ML_Operator *Amatrix,
				  int *aggr_count, int nvertices, 
				  int *aggr_index, int exp_Nrows, 
				  ML_Comm *comm, char *input_bdry,char *label);


/* ------------------------------------------------------------------------- */
/* set threshold for pruning matrix graph                                    */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_Threshold( ML_Aggregate *, double epsilon );
int ML_Aggregate_Reset_Threshold( ML_Aggregate * );

/* ------------------------------------------------------------------------- */
/* whether to smooth existing tentative prolongator                          */
/* ------------------------------------------------------------------------- */

extern int ML_Aggregate_Set_Flag_SmoothExistingTentativeP( ML_Aggregate *, int);
extern int ML_Aggregate_Get_Flag_SmoothExistingTentativeP( ML_Aggregate *);

/* ------------------------------------------------------------------------- */
/* damping factor for prolongator smoother                                   */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_DampingFactor( ML_Aggregate *, double factor );
int ML_Aggregate_Set_PSmootherType( ML_Aggregate *, int stype );

/* ------------------------------------------------------------------------- */
/* set up scheme to compute spectral radius of A at each level               */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_SpectralNormScheme_Calc( ML_Aggregate * );
int ML_Aggregate_Set_SpectralNormScheme_Anorm( ML_Aggregate * );

/* ------------------------------------------------------------------------- */
/* accessing aggregation information                                         */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Get_AggrCount( ML_Aggregate *, int level );
int ML_Aggregate_Get_AggrMap( ML_Aggregate *, int level, int**);
extern int ML_Gen_Blocks_Aggregates(ML_Aggregate *ag, int level, 
                                    int *nblocks, int **block_list);


/* ------------------------------------------------------------------------- */
/* set null space for the finest grid                                        */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_NullSpace(ML_Aggregate *, int, int, double *, int);
int ML_Aggregate_Scale_NullSpace(ML_Aggregate *ag, double *scale_vect,
				 int length);

int ML_Aggregate_Coarsen( ML_Aggregate *, ML_Operator *A,
                          ML_Operator **P, ML_Comm *comm );

/* ------------------------------------------------------------------------- */
/* functions for performing aggregation                                      */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_CoarsenMIS(ML_Aggregate *ml_ag,ML_Operator *Amatrix,
                            ML_Operator **Pmatrix, ML_Comm *comm);

int ML_Aggregate_CoarsenDomainDecomp(ML_Aggregate *ml_ag,
           ML_Operator *Amatrix, ML_Operator **Pmatrix, ML_Comm *comm);

int ML_Aggregate_CoarsenUncoupled(ML_Aggregate *ml_ag,
           ML_Operator *Amatrix, ML_Operator **Pmatrix, ML_Comm *comm);

int ML_Aggregate_CoarsenCoupled(ML_Aggregate *ml_ag,
           ML_Operator *Amatrix, ML_Operator **Pmatrix, ML_Comm *comm);

/* ------------------------------------------------------------------------- */
/* miscellaneous local functions                                             */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_LabelVertices(int, int *, char, char *,char *,int,
              int *, int *, int, int **, int, int **, int *, int *,
              int, int **, int *, int *, int **, int, ML_Comm *, int *);

int ML_Aggregate_UpdateVertexStates(int N_remaining_vertices,
              char vertex_state[], int recv_cnt, int recv_proc[], 
              int recv_leng[], int **recv_buf, int **recv_list, 
              int proc_flag[], int *NremainingRcvProcs, int send_cnt, 
              int send_proc[], int send_leng[], int **send_buf, 
              int *send_flag, USR_REQ *Request, ML_Comm *comm, int);

int ML_Aggregate_ExchangeBdry(double *data, void *comm);

int ML_Aggregate_ExchangeData(char *rbuf, char *sbuf, int N_neigh,
                              int *neigh, int *rleng, int *sleng,
                              int msgid, int datatype, ML_Comm *comm);

int ML_Aggregate_Print(ML_Aggregate *);
int ML_Aggregate_Print_Complexity(ML_Aggregate *);

void ML_CSR_MSR_ML_memorydata_Destroy(void *data);

#ifdef __cplusplus
}
#endif

#endif

