/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/

/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


/***********************************************************************

  Fortran wrappers so that Aztec may be called from fortran programs.

  ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "az_aztec.h"

extern void AZ_BROADCAST_F77(char *ptr, int *length, int proc_config[],int *action);
extern void AZ_CHECK_MSR_F77(int *bindx, int *N_update, int *N_external,
			  int option, int *proc_config);
extern void AZ_CHECK_VBR_F77(int *N_update, int *N_external, int *option, int bindx[],
			  int bnptr[], int cnptr[], int rnptr[], int proc_config[]);
extern void AZ_DEFAULTS_F77(int options[], double params[]);
extern void AZ_EXCHANGE_BDRY_F77(double x[], int data_org[], int proc_config[]);
extern void AZ_FIND_LOCAL_INDICES_F77(int *N_update, int bindx[], int update[],
			           int *external, int *N_external, int *mat_type,
                                   int bpntr[]);
extern void AZ_FIND_PROCS_FOR_EXTERNS_F77(int *N_update, int update[], int external[],
				       int *N_external, int proc_config[],
                                       int *extern_proc);
extern void AZ_FREE_MEMORY_F77(int *label);
extern void   AZ_GSUM_VEC_INT_F77(int vals[], int vals2[], int *length, int proc_config[]);
extern void   AZ_INIT_QUICK_FIND_F77(int list[], int *length, int *shift, int *bins);
extern void AZ_INVORDER_VEC_F77(double vector[], int data_org[], int update_index[],
			     int rpntr[], double newvec[]);
extern void   AZ_MATVEC_MULT_F77(double *val, int *indx, int *bindx, int *rpntr,
			      int *cpntr, int *bpntr, double *b, double *c,
                              int *exchange_flag, int *data_org);

extern void   AZ_MSR2VBR_F77(double val[], int indx[], int rnptr[], int cnptr[],
			  int bnptr[], int bindx[], int msr_bindx[], double msr_val[],
                          int *total_blk_rows, int *total_blk_cols, int *blk_space,
                          int *nz_space, int *blk_type);
extern void   AZ_ORDER_ELE_F77(int update_index[], int extern_index[], int *internal,
			    int *border, int *N_update, int msr_bindx[], int bindx[],
                            int extern_proc[], int *N_external, int *option,
                            int *m_type);
extern void   AZ_PRINT_ERROR_F77(int *error_code);
extern void   AZ_PROCESSOR_INFO_F77(int proc_config[]);
extern void   AZ_READ_MSR_MATRIX_F77(int update[], double *val, int *bindx, int *N_update,
				  int proc_config[]);
extern void   AZ_READ_UPDATE_F77(int *N_update_blks, int update_blks[], int proc_config[],
			      int *bigN, int *chunk, int *input_option);
extern void   AZ_REORDER_MATRIX_F77(int *N_update, int bindx[], double val[],
				 int update_index[], int extern_index[], int indx[],
                                 int rnptr[], int bnptr[], int *N_external,
                                 int cnptr[], int *option, int *mat_type);
extern void AZ_REORDER_VEC_F77(double vector[], int data_org[], int update_index[],
			    int rpntr[]);
extern void AZ_SET_COMM_F77(int proc_config[], MPI_AZComm comm);
extern void AZ_SET_PROC_CONFIG_F77(int proc_config[], MPI_AZComm *comm);
extern void   AZ_SET_MESSAGE_INFO_F77(int *N_external, int extern_index[],
			           int *N_update, int external[], int extern_proc[],
                            	   int update[], int update_index[], int proc_config[],
                                   int cnptr[], int *data_org, int *mat_type);
extern void AZ_SOLVE_F77(double x[], double b[], int options[], double params[],
		      int indx[], int bindx[], int rpntr[], int cpntr[], int bpntr[],
                      double val[], int data_org[], double status[], int proc_config[]);
extern void AZ_SORT_F77(int list[], int *N, int list2[], double list3[]);
extern void AZ_TRANSFORM_F77(int proc_config[],int external[],int bindx[], double val[],
			  int update[], int update_index[], int extern_index[],
                     int data_org[], int *N_update, int indx[], int bnptr[],
                     int rnptr[], int cnptr[], int *mat_type);
double AZ_GAVG_DOUBLE_F77(double *var, int proc_config[]);
double AZ_GDOT_F77(int *N, double r[], double z[],int proc_config[]);
double AZ_GMAX_DOUBLE_F77(double *var, int proc_config[]);
double AZ_GMAX_MATRIX_NORM_F77(double val[], int indx[], int bindx[], int rpntr[],
                            int cpntr[], int bpntr[], int proc_config[],
                            int data_org[]);
double AZ_GMAX_VEC_F77(int *N, double vec[], int proc_config[]);
double AZ_GMIN_DOUBLE_F77(double *var, int proc_config[]);
double AZ_GSUM_DOUBLE_F77(double *var, int proc_config[]);
double AZ_GVECTOR_NORM_F77(int *n, int *p, double *x, int *proc_config);

int AZ_CHECK_INPUT_F77(int data_org[], int options[], double params[],
                     int proc_config[]);
int AZ_FIND_INDEX_F77(int *key, int list[], int *length);
int AZ_GMAX_INT_F77(int *val, int proc_config[]);
int AZ_GMIN_INT_F77(int *val, int proc_config[]);
int AZ_GSUM_INT_F77(int *totals, int proc_config[]);
int AZ_QUICK_FIND_F77(int *key, int list[],int *length, int *shift, int bins[]);
extern MPI_AZComm *AZ_GET_COMM_F77(int proc_config[]);


int AZ_using_fortran = AZ_FALSE;

void AZ_BROADCAST_F77(char *ptr, int *length, int proc_config[],int *action)
{
  AZ_broadcast(ptr, *length, proc_config, *action);
}

int  AZ_CHECK_INPUT_F77(int data_org[], int options[], double params[],
                     int proc_config[])
{
  /* untouched */
  return(AZ_check_input(data_org, options, params, proc_config));
}

void AZ_CHECK_MSR_F77(int *bindx, int *N_update, int *N_external,
                   int option, int *proc_config)
{
  AZ_check_msr(bindx, *N_update, *N_external, option, proc_config);
}

void AZ_CHECK_VBR_F77(int *N_update, int *N_external, int *option, int bindx[],
                   int bnptr[], int cnptr[], int rnptr[], int proc_config[])
{
  AZ_check_vbr(*N_update, *N_external, *option, bindx, bnptr, cnptr, rnptr,
               proc_config);
}

void AZ_DEFAULTS_F77(int options[], double params[])
{
  /* untouched */
  AZ_defaults(options, params);
}

void AZ_EXCHANGE_BDRY_F77(double x[], int data_org[], int proc_config[])
{
  /* untouched */
  AZ_exchange_bdry(x, data_org, proc_config);
}

int  AZ_FIND_INDEX_F77(int *key, int list[], int *length)
{
  return(AZ_find_index(*key, list, *length));
}

void AZ_FIND_LOCAL_INDICES_F77(int *N_update, int bindx[], int update[],
                            int *external, int *N_external, int *mat_type,
                            int bpntr[])
{
  AZ_using_fortran = AZ_TRUE;
  AZ_find_local_indices(*N_update, bindx, update, &external,
                        N_external, *mat_type, bpntr);
  AZ_using_fortran = AZ_FALSE;
}

void AZ_FIND_PROCS_FOR_EXTERNS_F77(int *N_update, int update[], int external[],
                                int *N_external, int proc_config[],
                                int *extern_proc)
{

  AZ_using_fortran = AZ_TRUE;
  AZ_find_procs_for_externs(*N_update, update, external, *N_external,
                            proc_config, &extern_proc);
  AZ_using_fortran = AZ_FALSE;
}

void AZ_FREE_MEMORY_F77(int *label)
{
  AZ_free_memory(*label);
}

double AZ_GAVG_DOUBLE_F77(double *var, int proc_config[])
{
  return(AZ_gavg_double(*var, proc_config));
}

double AZ_GDOT_F77(int *N, double r[], double z[],int proc_config[])
{
  return(AZ_gdot(*N, r, z,proc_config));
}

double AZ_GMAX_DOUBLE_F77(double *var, int proc_config[])
{
  return(AZ_gmax_double(*var, proc_config));
}

int AZ_GMAX_INT_F77(int *val, int proc_config[])
{
  return(AZ_gmax_int(*val, proc_config));
}

double AZ_GMAX_MATRIX_NORM_F77(double val[], int indx[], int bindx[], int rpntr[],
                            int cpntr[], int bpntr[], int proc_config[],
                            int data_org[])
{
  /* unchanged */
  return(AZ_gmax_matrix_norm(val, indx, bindx, rpntr, cpntr, bpntr,
                             proc_config, data_org));
}

double AZ_GMAX_VEC_F77(int *N, double vec[], int proc_config[])
{
  return(AZ_gmax_vec(*N, vec, proc_config));
}

double AZ_GMIN_DOUBLE_F77(double *var, int proc_config[])
{
  return(AZ_gmin_double(*var, proc_config));
}

int AZ_GMIN_INT_F77(int *val, int proc_config[])
{
  return(AZ_gmin_int(*val, proc_config));
}

double AZ_GSUM_DOUBLE_F77(double *var, int proc_config[])
{
  return(AZ_gsum_double(*var, proc_config));
}

int    AZ_GSUM_INT_F77(int *totals, int proc_config[])
{
  return(AZ_gsum_int(*totals, proc_config));
}

void   AZ_GSUM_VEC_INT_F77(int vals[], int vals2[], int *length, int proc_config[])
{
  AZ_gsum_vec_int(vals, vals2, *length, proc_config);
}

double AZ_GVECTOR_NORM_F77(int *n, int *p, double *x, int *proc_config)
{
  return(AZ_gvector_norm(*n, *p, x, proc_config));
}

void   AZ_INIT_QUICK_FIND_F77(int list[], int *length, int *shift, int *bins)
{
  AZ_init_quick_find(list, *length, shift, bins);
}

void AZ_INVORDER_VEC_F77(double vector[], int data_org[], int update_index[],
	int rpntr[], double newvec[])
{
  AZ_invorder_vec(vector, data_org, update_index, rpntr, newvec);
}

void   AZ_MATVEC_MULT_F77(double *val, int *indx, int *bindx, int *rpntr,
                       int *cpntr, int *bpntr, double *b, double *c,
                       int *exchange_flag, int *data_org)
{
  AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, b, c,
                 *exchange_flag, data_org);
}

void   AZ_MSR2VBR_F77(double val[], int indx[], int rnptr[], int cnptr[],
                   int bnptr[], int bindx[], int msr_bindx[], double msr_val[],
                   int *total_blk_rows, int *total_blk_cols, int *blk_space,
                   int *nz_space, int *blk_type)
{
  AZ_msr2vbr(val, indx, rnptr, cnptr, bnptr, bindx, msr_bindx, msr_val,
             *total_blk_rows, *total_blk_cols, *blk_space, *nz_space,
             *blk_type);
}

void   AZ_ORDER_ELE_F77(int update_index[], int extern_index[], int *internal,
                     int *border, int *N_update, int msr_bindx[], int bindx[],
                     int extern_proc[], int *N_external, int *option,
                     int *m_type)
{
  AZ_order_ele(update_index, extern_index, internal, border, *N_update,
               msr_bindx, bindx, extern_proc, *N_external, *option, *m_type);
}

void   AZ_PRINT_ERROR_F77(int *error_code)
{
  AZ_print_error(*error_code);
}

void   AZ_PROCESSOR_INFO_F77(int proc_config[])
{
  /* no change */
  AZ_processor_info(proc_config);
}

int    AZ_QUICK_FIND_F77(int *key, int list[],int *length, int *shift, int bins[])
{
  return(AZ_quick_find(*key, list,*length, *shift, bins));
}

void   AZ_READ_MSR_MATRIX_F77(int update[], double *val, int *bindx, int *N_update,
                           int proc_config[])
{
  AZ_using_fortran = AZ_TRUE;
  AZ_read_msr_matrix(update, &val, &bindx, *N_update, proc_config);
  AZ_using_fortran = AZ_FALSE;
}

void   AZ_READ_UPDATE_F77(int *N_update_blks, int update_blks[], int proc_config[],
                       int *bigN, int *chunk, int *input_option)
{
  AZ_using_fortran = AZ_TRUE;
  AZ_read_update(N_update_blks, &update_blks, proc_config, *bigN, *chunk,
                 *input_option);
  AZ_using_fortran = AZ_FALSE;
}

void   AZ_REORDER_MATRIX_F77(int *N_update, int bindx[], double val[],
                          int update_index[], int extern_index[], int indx[],
                          int rnptr[], int bnptr[], int *N_external,
                          int cnptr[], int *option, int *mat_type)
{
  AZ_reorder_matrix(*N_update, bindx, val, update_index, extern_index, indx,
		    rnptr, bnptr, *N_external, cnptr, *option, *mat_type);
}

void AZ_REORDER_VEC_F77(double vector[], int data_org[], int update_index[],
        int rpntr[])
{
  AZ_reorder_vec(vector, data_org, update_index, rpntr);
}

MPI_AZComm *AZ_GET_COMM_F77(int proc_config[]) { return(AZ_get_comm(proc_config));}
void AZ_SET_COMM_F77(int proc_config[], MPI_AZComm comm) {
   AZ_set_comm(proc_config, comm); }
void AZ_SET_PROC_CONFIG_F77(int proc_config[], MPI_AZComm *comm) {
   AZ_set_proc_config(proc_config, *comm); }

void   AZ_SET_MESSAGE_INFO_F77(int *N_external, int extern_index[],
                            int *N_update, int external[], int extern_proc[],
                            int update[], int update_index[], int proc_config[],
                            int cnptr[], int *data_org, int *mat_type)
{
  AZ_using_fortran = AZ_TRUE;
  AZ_set_message_info(*N_external, extern_index,
                      *N_update, external,
                      extern_proc, update,
                      update_index, proc_config,
                      cnptr, &data_org, *mat_type);
  AZ_using_fortran = AZ_FALSE;
}

void AZ_SOLVE_F77(double x[], double b[], int options[], double params[],
               int indx[], int bindx[], int rpntr[], int cpntr[], int bpntr[],
               double val[], int data_org[], double status[], int proc_config[])
{
  /* no change */
  AZ_solve(x, b, options, params, indx, bindx, rpntr, cpntr, bpntr,
           val, data_org, status, proc_config);
}

void   AZ_SORT_F77(int list[], int *N, int list2[], double list3[])
{
  AZ_sort(list, *N, list2, list3); 
}

void   AZ_TRANSFORM_F77(int proc_config[],int external[],int bindx[], double val[],
                     int update[], int update_index[], int extern_index[],
                     int data_org[], int *N_update, int indx[], int bnptr[],
                     int rnptr[], int cnptr[], int *mat_type)
{
  AZ_using_fortran = AZ_TRUE;
  AZ_transform(proc_config,&external,bindx, val, update, &update_index,
               &extern_index, &data_org, *N_update, indx, bnptr, rnptr, &cnptr,
               *mat_type);
  AZ_using_fortran = AZ_FALSE;
}

