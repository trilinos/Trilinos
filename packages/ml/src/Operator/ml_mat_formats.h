/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of matrix-format specific stuff                          */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Ray Tuminaro (SNL)           */
/* Date          : April,    1999                                       */
/* ******************************************************************** */

#ifndef _MLMATFORMATS_
#define _MLMATFORMATS_


/* ******************************************************************** */
/* Structure used for ML_MSR_getrows and ML_CSR_getrows                 */
/*                                                                      */
/* For CSR matrices, we have                                            */
/*    A(i,j) where 0 <= i < N |    j = columns[k] , A(i,j) = values[k]  */
/*           and   0 <= j < N |    where  rowptr[i] <= k <rowptr[i+1]   */
/*                                                                      */
/* For MSR matrices, we have                                            */
/*    A(i,j) where 0 <= i < N,|    j = columns[k] , A(i,j) = values[k]  */
/*                 i != j ,   |    where  columns[i] <= k < columns[i+1]*/
/*             and 0 <= j < N |                                         */
/*    A(i,i) where 0 <= i < N,|    A(i,j) = values[i]                   */
/* -------------------------------------------------------------------- */

struct ML_CSR_MSRdata 
{
   int    *columns, *rowptr;
   double *values;
                              /*************************************/
  int    Nnz, Nrows, Ncols;   /* Haim's addition. Convenient for GGB */
                              /* implementation of prolongator.      */
                              /* These fields are not normally filled */
                              /* in the rest of ml. */
};

struct ML_vbrdata 
{
   int    *bindx, *bpntr, *cpntr, *rpntr, *indx;
   double *val;
};
#include "ml_common.h"
#include "ml_comminfoop.h"

typedef struct ML_Matrix_DCSR_Struct
{
   int           ML_id;
   int           mat_n;
   int           *mat_ia;
   int           *mat_ja;
   double        *mat_a;
   ML_Comm       *comm;
   ML_CommInfoOP *comminfo;

} ML_Matrix_DCSR;

#include "ml_operator.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


extern void ML_RECUR_CSR_MSRdata_Destroy(ML_Operator *matrix);
extern void ML_OnlyFreeTopLevelDataPtr(void *data);
extern void ML_CSR_MSRdata_Destroy(void *data);
extern void ML_CSR_MSRdata_Destroy_StructOnly(void *data);

extern void ML_restricted_MSR_mult(ML_Operator *matrix, int Nrows,
                                   double b[], double c[], int Nsend);

extern void ML_Scale_CSR(ML_Operator *input_matrix,
                         double scale_factors[], int mult_or_divide);

extern int CSR_getrows(void *data,int N_requested_rows,int requested_rows[],
                       int allocated_space, int columns[], double values[],
                       int row_lengths[]);
extern int sCSR_getrows(void *data,int N_requested_rows,int requested_rows[],
                       int allocated_space, int columns[], double values[],
                       int row_lengths[]);
extern int cCSR_getrows(void *data, int N_requested_rows, int requested_rows[],
                       int allocated_space, int columns[], double values[],
                       int row_lengths[]);
extern int cCSR_trans_matvec(void *Amat_in, int ilen, double p[], int olen, 
			     double ap[]);
extern int cCSR_matvec(void *Amat_in, int ilen, double p[], int olen, 
		       double ap[]);

extern int CSR_get_ones_rows(void *data, int N_requested_rows, 
			     int requested_rows[], int allocated_space, 
			     int columns[], double values[], int row_lengths[]);


extern int MSR_get_ones_rows(void *data, int N_requested_rows, int requested_rows[],
                       int allocated_space, int columns[], double values[],
                       int row_lengths[]);
extern int MSR_getrows(void *data, int N_requested_rows, int requested_rows[],
                       int allocated_space, int columns[], double values[],
                       int row_lengths[]);

extern int MSR_matvec(void *Amat, int, double p[], int, double ap[]);

#ifdef WKC
/* WKC -- double * happen to be Epetra_MultiVectors in cognito */
extern int MSR_matvec_WKC(void *Amat, int, double *p, int, double *ap);
#endif

extern int CSR_densematvec(void *Amat_in, int ilen, double p[], int olen, double ap[]);
extern int CSR_matvec(void *Amat, int, double p[], int, double ap[]);
extern int sCSR_trans_matvec(void *Amat, int, double p[], int, double ap[]);
extern int sCSR_matvec(void *Amat, int, double p[], int, double ap[]);
extern int CSR_ones_matvec(void *Amat, int, double p[], int, double ap[]);
extern int localCSR_matvec(void *Amat_in, int ilen, double p[], int olen, 
                           double ap[]);


extern int VBR_cnst_blk_getrows(void *data, int N_requested_rows,
                                int requested_rows[], int allocated_space,
                                int columns[], double values[],
                                int row_lengths[]);

extern int VECTOR_getrows(void *data,int N_requested_rows,int requested_rows[],
                          int allocated_space, int columns[], double values[],
                          int row_lengths[]);
extern int ML_MSR2CSR(struct ML_CSR_MSRdata *csr_data, int Nrows,
                          int *Ncolumns);  

extern int  ML_Matrix_DCSR_Create( ML_Matrix_DCSR ** );
extern void ML_Matrix_DCSR_Destroy( ML_Matrix_DCSR ** );
extern int  ML_Matrix_DCSR_Set( ML_Matrix_DCSR *,int,int*,int*,double*);
extern int  ML_Matrix_DCSR_Set_Comm(ML_Matrix_DCSR*,ML_CommInfoOP*,ML_Comm*);
extern int  ML_Matrix_DCSR_Getrow(void*,int,int*,int,int*,double*,int*);
extern int  ML_Matrix_DCSR_Matvec(void*,int,double*,int,double*);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif
#endif

