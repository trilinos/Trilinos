/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* data structure to hold AMG information                                    */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL)                                       */
/* Date          : October, 2000                                             */
/* ************************************************************************* */

#ifndef __MLAMGH__
#define __MLAMGH__

#include <stdio.h>
/* #include <stdlib.h> */
#include "ml_common.h"
#include "ml_defs.h"

/* ************************************************************************* */
/* local defines                                                             */
/* ------------------------------------------------------------------------- */

#define ML_AMG_MIS              21
#define ML_AMG_SCALAR            0
#define ML_AMG_SYSTEM_UNKNOWN    1
#define ML_AMG_SYSTEM_NODAL      2
#define ML_AMG_SM_JACOBI        11
#define ML_AMG_SM_GS            12
#define ML_AMG_SM_SGS           13
#define ML_AMG_SM_VBJACOBI      14
#define ML_AMG_SM_VBGS          15
#define ML_AMG_SM_VBSGS         16
#define ML_AMG_SM_ASCHWARZ      17
#define ML_AMG_SM_MSCHWARZ      18
#define ML_AMG_SM_SUPERLU       19

/* ************************************************************************* */
/* other ML include files                                                    */
/* ------------------------------------------------------------------------- */

#include "ml_defs.h"
#include "ml_comm.h"
#include "ml_memory.h"
#include "ml_operator.h"

/* ************************************************************************* */
/* definition of the AMG structure                                           */
/* ------------------------------------------------------------------------- */

typedef struct ML_AMG_Struct
{
   int    ML_id;
   double print_flag;
   int    max_coarse_size;             /** maximum size of coarsest grid */
   double threshold;                   /** for pruning matrix            */ 
   double curr_threshold;              /** adjusted for levels           */
   int    coarsen_scheme;              /** MIS                           */
   int    amg_scheme;                  /** scalar(0),unknown(1),system(2)*/
   int    num_PDE_eqns;                /** block size                    */
   int    *blk_info;                   /** store dof information         */
   int    max_levels;                  /** maximum number of levels      */
   int    begin_level;                 /** finest grid level             */
   int    cur_level;                   /** temporary variable            */
   double operator_complexity;         /** sum of nnz for all A's        */
   double fine_complexity;             /** nnz of the finest A           */

   ML_Operator *fine_Amat;
   int    presmoother_type;
   int    postsmoother_type;
   int    presmoother_ntimes;
   int    postsmoother_ntimes;
   double presmoother_jacobiwt;
   double postsmoother_jacobiwt;
   int    coarse_solver_type;
   int    coarse_solver_ntimes;
   double coarse_solver_jacobiwt;
   int    *pre_aztec_options;
   double *pre_aztec_params;
   int    *pre_aztec_proc_config;
   double *pre_aztec_status;
   void   (*pre_function)(int);
  /* Trying to avoid including az_aztec.h
   void   (*pre_function)(double *, int *, int *, double *,
                          struct AZ_MATRIX_STRUCT  *,
                          struct AZ_PREC_STRUCT *);
  */
   int    *post_aztec_options;
   double *post_aztec_params;
   int    *post_aztec_proc_config;
   double *post_aztec_status;
   void   (*post_function)(int);
  /* Trying to avoid including az_aztec.h
   void   (*post_function)(double *, int *, int *, double *,
                           struct AZ_MATRIX_STRUCT  *,
                           struct AZ_PREC_STRUCT *);
  */
} ML_AMG;

/* ************************************************************************* */
/* ************************************************************************* */
/* functions to manipulate the AMG data structure                            */
/* ------------------------------------------------------------------------- */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" 
{
#endif
#endif

/* ------------------------------------------------------------------------- */
/* constructor/destructor and level control                                  */
/* ------------------------------------------------------------------------- */

extern int  ML_AMG_Create( ML_AMG ** );
extern int  ML_AMG_Destroy( ML_AMG ** );

extern int  ML_AMG_Set_OutputLevel( ML_AMG *amg, double outputlevel );

extern int  ML_AMG_Set_MaxLevels( ML_AMG *amg, int level );
extern int  ML_AMG_Set_CurrentLevel( ML_AMG *amg, int level );
extern int  ML_AMG_Set_StartLevel( ML_AMG *amg, int level );

/* ------------------------------------------------------------------------- */
/* control when to stop coarsening                                           */
/* ------------------------------------------------------------------------- */

extern int  ML_AMG_Set_MaxCoarseSize( ML_AMG *amg, int size );

/* ------------------------------------------------------------------------- */
/* different AMG scheme (scalar, unknown, system)                            */
/* ------------------------------------------------------------------------- */

extern int  ML_AMG_Set_AMGScheme_Scalar( ML_AMG *amg  );
extern int  ML_AMG_Set_AMGScheme_SystemUnknown( ML_AMG *amg, int  );
extern int  ML_AMG_Set_AMGScheme_SystemNodal( ML_AMG *amg, int  );

/* ------------------------------------------------------------------------- */
/* different parallel coarsening schemes                                     */
/* ------------------------------------------------------------------------- */

extern int  ML_AMG_Set_CoarsenScheme_MIS( ML_AMG *amg  );

/* ------------------------------------------------------------------------- */
/* set threshold for pruning matrix graph                                    */
/* ------------------------------------------------------------------------- */

extern int  ML_AMG_Set_Threshold( ML_AMG *amg, double epsilon );

/* ------------------------------------------------------------------------- */
/* functions for performing coarsening                                       */
/* ------------------------------------------------------------------------- */

extern int  ML_AMG_Coarsen(ML_AMG *amg, ML_Operator *Amatrix,
                           ML_Operator **Pmatrix, ML_Comm *comm);
extern int  ML_AMG_CoarsenMIS(ML_AMG *ml_amg,ML_Operator *Amatrix,
                              ML_Operator **Pmatrix, ML_Comm *comm);

/* ------------------------------------------------------------------------- */
/* miscellaneous functions                                                   */
/* ------------------------------------------------------------------------- */

extern int  ML_AMG_Set_Smoother(ML_AMG *ml_amg,int smoother_type,int pre_or_post,
                        ML_Operator *Amatrix,int ntimes, double weight);

extern int  ML_AMG_Set_SmootherAztec(ML_AMG *ml_amg, int pre_or_post,
                       int *options, double *params, int *proc_config,
				     double *status, void (*post_function)(int)
				     /* Trying to avoid including az_aztec.h 
                       void (*post_function)(double *,int *,int *,double *,
		       struct AZ_MATRIX_STRUCT*,struct AZ_PREC_STRUCT*)*/);

extern int  ML_AMG_Set_CoarseSolve(ML_AMG *ml_amg, int solve_type, int ntimes, 
                                   double weight);

/* ------------------------------------------------------------------------- */
/* miscellaneous functions                                                   */
/* ------------------------------------------------------------------------- */

extern int  ML_AMG_Print(ML_AMG *);
extern int  ML_AMG_Print_Complexity(ML_AMG *);

extern int  ML_AMG_LabelVertices(int nvertices, int *vlist, char,
                  char *state, char *vtype, int nvert, int *rowptr, 
                  int *columns, int mypid, int **proclist, int Nneigh,
                  int **sndbuf, int *neigh, int *sndleng, int Nneigh2,
                  int **rcvbuf, int *neigh2, int *rcvleng, int **recvlist, 
                  ML_Comm *comm, int *CF_array);


int ML_AMG_GetCommInfo(ML_CommInfoOP *mat_comm, int Nrows, int *A_Nneigh, 
           int **A_neigh, int ***A_sendlist, int ***A_recvlist, 
           int ***A_sndbuf, int ***A_rcvbuf, int **A_sndleng, 
           int **A_rcvleng, int *Nghost);

int ML_AMG_UpdateVertexStates(int N_remaining_vertices, char vertex_state[], 
           int recv_cnt, int recv_proc[], int recv_leng[], int **recv_buf, 
           int **recv_list, int proc_flag[], int *NremainingRcvProcs, 
           int send_cnt, int send_proc[], int send_leng[], int **send_buf, 
           int *send_flag, USR_REQ *Request, ML_Comm *comm, int msgtype);


int ML_AMG_CompatibleRelaxation(int *CF_array,
           ML_Operator *Amat, int *Ncoarse, int limit);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif

