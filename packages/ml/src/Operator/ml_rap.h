/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML RAP data structure                             */
/* ******************************************************************** */
/* Author        : Raymond Tuminaro (SNL)                               */
/* Date          : February, 1999                                       */
/* ******************************************************************** */

#ifndef _MLRAP_
#define _MLRAP_

#include "ml_defs.h"
#include "ml_struct.h"
#include "ml_comminfoagx.h"

/* ******************************************************************** */
/* external function proto-types                                        */
/* -------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

extern void ML_convert_data_org(ML_Operator *mat, int d_org[], int *rcv,
                                int *remap, int leng, int add_or_not);

extern void ML_add_appended_rows(ML_CommInfoOP *comm_info, 
                                 ML_Operator *matrix, int orig_rows, 
                                 int total_rcvd, int appended_nzs);

extern void ML_back_to_local(ML_Operator *imatrix, ML_Operator *omatrix, 
                             int max_per_proc);

extern void ML_back_to_csrlocal(ML_Operator *imatrix, ML_Operator *omatrix,
        int max_per_proccomm);

extern void ML_exchange_rows(ML_Operator *orig, ML_Operator **appended, 
                             ML_CommInfoOP *comm_info);

extern void ML_expand_accum(int accum_size, int **accum_col, 
                            double **accum_val, int Ncols);

extern void ML_get_matrix_row(ML_Operator *input_matrix,int N_requested_rows,
                              int requested_rows[], int *allocated_space, 
                              int **columns, double **values, 
                              int row_lengths[], int index);

extern void ML_matmat_mult(ML_Operator *Amat, ML_Operator *Bmat, 
                           ML_Operator **Cmat);

extern void ML_getrow_matvec(ML_Operator *matrix, double *vec, 
                             int Nvec, double *ovec, int *Novec);

extern void ML_rap(ML_Operator *R2mat, ML_Operator *A2mat, 
                   ML_Operator *P2mat, ML_Operator *Result);

extern void ML_rap_check(ML *ml, ML_Operator *RAP, ML_Operator *R,
                         ML_Operator *A, ML_Operator *P, int iNvec,
                         int oNvec);

extern void ML_set_message_info(int N_external, int external[], 
                                int max_per_proc, ML_Operator *omatrix);

extern void ML_sum_duplicates(int accum_col[],double accum_val[],int *Ncols);

#ifdef __cplusplus
}
#endif

#define ML_allocate(i)    malloc((i))
#define ML_free(i)        { free(i); i = NULL; }

#define ML_matrix_type         0
#define ML_N_internal          1
#define ML_N_border            2
#define ML_N_external          3
#define ML_N_int_blk           4
#define ML_N_bord_blk          5
#define ML_N_ext_blk           6
#define ML_N_neigh             7
#define ML_total_send          8
#define ML_name                9
#define ML_internal_use        10
#define ML_N_rows              11
#define ML_neighbors           12
#define ML_rec_length          (12 +   ML_MAX_NEIGHBORS)
#define ML_send_length         (12 + 2*ML_MAX_NEIGHBORS)
#define ML_send_list           (12 + 3*ML_MAX_NEIGHBORS)
#define ML_FALSE               0
#define ML_TRUE                1

#ifdef __cplusplus
extern "C"
{
#endif

extern double ML_gsum_double(double val, ML_Comm *comm);

extern void ML_az_sort(int list[], int N, int list2[], double list3[]);

extern void ML_gsum_vec_int(int vals[],int vals2[],int length,ML_Comm *comm);

extern void ML_rm_duplicates(int array[], int *N);

extern void ML_splitup_big_msg(int num_neighbors,char *ibuffer,char *obuffer,
                               unsigned int element_size,int *start_send_proc,
                               int *actual_send_length,int *actual_recv_length,
                               int *proc_num_neighbor, int type, 
                               int *total_num_recv, ML_Comm *comm);
extern double ML_gdot(int N, double r[], double z[], ML_Comm *comm);

extern int ML_gmax_int(int val, ML_Comm *comm);

#ifdef __cplusplus
}
#endif

#define max(x,y) (( x > y ) ? x : y)     /* max function  */

#ifdef __cplusplus
extern "C"
{
#endif

extern double ddot_(int *n1, double *v1, int *dum11, double *v2, int *dum21);

extern int    ML_find_index(int key, int list[], int length);

#ifdef __cplusplus
}
#endif

#endif
