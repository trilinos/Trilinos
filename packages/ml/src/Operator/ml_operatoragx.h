/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Data structure to hold the grid transfer operator.                   */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : April, 1998                                          */
/* ******************************************************************** */

#ifndef _MLOPERATORAGX_
#define _MLOPERATORAGX_

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/* ******************************************************************** */
/* these include are for (1) memory allocation, (2) communication       */
/* buffer access, and (3) communication (recv, send)                    */
/* -------------------------------------------------------------------- */

#include "ml_common.h"
#include "ml_memory.h"
#include "ml_comm.h"
#include "ml_comminfoagx.h"
#include "ml_struct.h"

/* ******************************************************************** */
/*  Nlocal_rows, local_ia, local_ja, local_a :                          */
/*        variables and arrays (in CSR format) to store local operator  */
/*  Nremote_rows, remote_ia, remote_ja, remote_a :                      */
/*        variables and arrays (in CSR format) to store remote operator */
/*  restrict_wgts : normalization factors for restriction               */
/*  ext arrays    : for temporary storage of interprocessor data        */
/*  node_flag     : for registering already processed fine nodes        */
/*  com           : processor communication pattern                     */
/* -------------------------------------------------------------------- */

typedef struct ML_OperatorAGX_Struct
{
   int    ML_id;
   int    proc_id, num_procs, step;
   int    Nlocal_rows, *local_ia, *local_ja;
   int    Nremote_rows, *remote_ia, *remote_ja;
   double *local_a, *remote_a, *restrict_wgts, *remote_restrict_wgts;
   int    ext_cnt, *ext_ia, *ext_ja;
   int    ext2_cnt, *ext2_ptr, *ext2_index;
   double *ext_a, *ext2_a;
   int    *fnode_flag, fnode_leng;
   int    *coarse_bdry_list, coarse_bdry_leng;
   int    AGX_stride;
   ML_Comm *AGX_comm;
   ML_CommInfoAGX *com;

} ML_OperatorAGX;


/* ******************************************************************** */
/* functions to manipulate the Operator structure                       */
/* -------------------------------------------------------------------- */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif

extern int ML_OperatorAGX_Clean_Getrows(ML_Operator **);
extern int ML_OperatorAGX_Create(ML_OperatorAGX **);
extern void ML_Operator2AGX_Destroy(void *);
extern int ML_OperatorAGX_Destroy(ML_OperatorAGX **);
extern int ML_OperatorAGX_Print(ML_OperatorAGX *);
extern int ML_OperatorAGX_Restrict(ML_Operator *,int,double *,int,double*);
extern int ML_OperatorAGX_Prolongate(ML_Operator *,int,double*,int,double*);
extern int ML_OperatorAGX_Getrows(ML_Operator *data, int N_requested_rows,
              int requested_rows[], int allocated_space, int columns[],
              double values[], int row_lengths[]);
extern int ML_OperatorAGX_Getcols(ML_Operator *data, int N_requested_rows,
              int requested_rows[], int allocated_space, int columns[],
	      double values[], int row_lengths[]);
extern int ML_OperatorAGX_Gen_ComminfoOp(ML_OperatorAGX *vop, ML_Operator *Rmat,
        ML_Operator *Pmat);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
