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

extern void az_broadcast_(char *ptr, int *length, int proc_config[],int *action);
extern void az_check_msr_(int *bindx, int *N_update, int *N_external,
			  int option, int *proc_config);
extern void az_check_vbr_(int *N_update, int *N_external, int *option, int bindx[],
			  int bnptr[], int cnptr[], int rnptr[], int proc_config[]);
extern void az_defaults_(int options[], double params[]);
extern void az_exchange_bdry_(double x[], int data_org[], int proc_config[]);
extern void az_find_local_indices_(int *N_update, int bindx[], int update[],
			           int *external, int *N_external, int *mat_type,
                                   int bpntr[]);
extern void az_find_procs_for_externs_(int *N_update, int update[], int external[],
				       int *N_external, int proc_config[],
                                       int *extern_proc);
extern void az_free_memory_(int *label);
extern void   az_gsum_vec_int_(int vals[], int vals2[], int *length, int proc_config[]);
extern void   az_init_quick_find_(int list[], int *length, int *shift, int *bins);
extern void az_invorder_vec_(double vector[], int data_org[], int update_index[],
			     int rpntr[], double newvec[]);
extern void   az_matvec_mult_(double *val, int *indx, int *bindx, int *rpntr,
			      int *cpntr, int *bpntr, double *b, double *c,
                              int *exchange_flag, int *data_org);

extern void   az_msr2vbr_(double val[], int indx[], int rnptr[], int cnptr[],
			  int bnptr[], int bindx[], int msr_bindx[], double msr_val[],
                          int *total_blk_rows, int *total_blk_cols, int *blk_space,
                          int *nz_space, int *blk_type);
extern void   az_order_ele_(int update_index[], int extern_index[], int *internal,
			    int *border, int *N_update, int msr_bindx[], int bindx[],
                            int extern_proc[], int *N_external, int *option,
                            int *m_type);
extern void   az_print_error_(int *error_code);
extern void   az_processor_info_(int proc_config[]);
extern void   az_read_msr_matrix_(int update[], double *val, int *bindx, int *N_update,
				  int proc_config[]);
extern void   az_read_update_(int *N_update_blks, int update_blks[], int proc_config[],
			      int *bigN, int *chunk, int *input_option);
extern void   az_reorder_matrix_(int *N_update, int bindx[], double val[],
				 int update_index[], int extern_index[], int indx[],
                                 int rnptr[], int bnptr[], int *N_external,
                                 int cnptr[], int *option, int *mat_type);
extern void az_reorder_vec_(double vector[], int data_org[], int update_index[],
			    int rpntr[]);
extern void az_set_comm_(int proc_config[], MPI_AZComm comm);
extern void az_set_proc_config_(int proc_config[], MPI_AZComm *comm);
extern void   az_set_message_info_(int *N_external, int extern_index[],
			           int *N_update, int external[], int extern_proc[],
                            	   int update[], int update_index[], int proc_config[],
                                   int cnptr[], int *data_org, int *mat_type);
extern void az_solve_(double x[], double b[], int options[], double params[],
		      int indx[], int bindx[], int rpntr[], int cpntr[], int bpntr[],
                      double val[], int data_org[], double status[], int proc_config[]);
extern void az_sort_(int list[], int *N, int list2[], double list3[]);
extern void az_transform_(int proc_config[],int external[],int bindx[], double val[],
			  int update[], int update_index[], int extern_index[],
                     int data_org[], int *N_update, int indx[], int bnptr[],
                     int rnptr[], int cnptr[], int *mat_type);
extern void az_broadcast__(char *ptr, int *length, int proc_config[],int *action);
extern void az_check_msr__(int *bindx, int *N_update, int *N_external,
			   int option, int *proc_config);
extern void az_check_vbr__(int *N_update, int *N_external, int *option, int bindx[],
			   int bnptr[], int cnptr[], int rnptr[], int proc_config[]);
extern void az_defaults__(int options[], double params[]);
extern void az_exchange_bdry__(double x[], int data_org[], int proc_config[]);
extern void az_find_local_indices__(int *N_update, int bindx[], int update[],
				    int *external, int *N_external, int *mat_type,
                            	    int bpntr[]);
extern void az_find_procs_for_externs__(int *N_update, int update[], int external[],
				        int *N_external, int proc_config[],
                                	int *extern_proc);
extern void az_free_memory__(int *label);
extern void   az_gsum_vec_int__(int vals[], int vals2[],int *length, int proc_config[]);
extern void   az_init_quick_find__(int list[], int *length, int *shift, int *bins);
extern void az_invorder_vec__(double vector[], int data_org[], int update_index[],
			      int rpntr[], double newvec[]);
extern void   az_matvec_mult__(double *val, int *indx, int *bindx, int *rpntr,
			       int *cpntr, int *bpntr, double *b, double *c,
                       	       int *exchange_flag, int *data_org);
extern void   az_msr2vbr__(double val[], int indx[], int rnptr[], int cnptr[],
			   int bnptr[], int bindx[], int msr_bindx[], double msr_val[],
                   	   int *total_blk_rows, int *total_blk_cols, int *blk_space,
                   	   int *nz_space, int *blk_type);
extern void   az_order_ele__(int update_index[], int extern_index[], int *internal,
			     int *border, int *N_update, int msr_bindx[], int bindx[],
                     	     int extern_proc[], int *N_external, int *option,
                             int *m_type);
extern void   az_print_error__(int *error_code);
extern void   az_processor_info__(int proc_config[]);
extern void   az_read_msr_matrix__(int update[],double *val, int *bindx, int *N_update,
				   int proc_config[]);
extern void   az_read_update__(int *N_update_blks,int update_blks[], int proc_config[],
			       int *bigN, int *chunk, int *input_option);
extern void   az_reorder_matrix__(int *N_update, int bindx[], double val[],
				  int update_index[], int extern_index[], int indx[],
                          	  int rnptr[], int bnptr[], int *N_external,
                          	  int cnptr[], int *option, int *mat_type);
extern void az_reorder_vec__(double vector[], int data_org[], int update_index[],
			     int rpntr[]);
extern void az_set_comm__(int proc_config[], MPI_AZComm comm);
extern void az_set_proc_config__(int proc_config[], MPI_AZComm *comm);
extern void   az_set_message_info__(int *N_external, int extern_index[],
				    int *N_update, int external[], int extern_proc[],
                            	    int update[], int update_index[], int proc_config[],
                            	    int cnptr[], int *data_org, int *mat_type);
extern void az_solve__(double x[], double b[], int options[], double params[],
               int indx[], int bindx[], int rpntr[], int cpntr[], int bpntr[],
               double val[], int data_org[], double status[], int proc_config[]);
extern void   az_sort__(int list[], int *N, int list2[], double list3[]);
extern void   az_transform__(int proc_config[],int external[],int bindx[], double val[],
			     int update[], int update_index[], int extern_index[],
                             int data_org[], int *N_update, int indx[], int bnptr[],
                             int rnptr[], int cnptr[], int *mat_type);
double az_gavg_double_(double *var, int proc_config[]);
double az_gdot_(int *N, double r[], double z[],int proc_config[]);
double az_gmax_double_(double *var, int proc_config[]);
double az_gmax_matrix_norm_(double val[], int indx[], int bindx[], int rpntr[],
                            int cpntr[], int bpntr[], int proc_config[],
                            int data_org[]);
double az_gmax_vec_(int *N, double vec[], int proc_config[]);
double az_gmin_double_(double *var, int proc_config[]);
double az_gsum_double_(double *var, int proc_config[]);
double az_gvector_norm_(int *n, int *p, double *x, int *proc_config);
double az_gavg_double__(double *var, int proc_config[]);
double az_gdot__(int *N, double r[], double z[],int proc_config[]);
double az_gmax_double__(double *var, int proc_config[]);
double az_gmax_matrix_norm__(double val[], int indx[], int bindx[], int rpntr[],
                            int cpntr[], int bpntr[], int proc_config[],
                            int data_org[]);
double az_gmax_vec__(int *N, double vec[], int proc_config[]);
double az_gmin_double__(double *var, int proc_config[]);
double az_gsum_double__(double *var, int proc_config[]);
double az_gvector_norm__(int *n, int *p, double *x, int *proc_config);

int az_check_input_(int data_org[], int options[], double params[],
                     int proc_config[]);
int az_find_index_(int *key, int list[], int *length);
int az_gmax_int_(int *val, int proc_config[]);
int az_gmin_int_(int *val, int proc_config[]);
int az_gsum_int_(int *totals, int proc_config[]);
int az_quick_find_(int *key, int list[],int *length, int *shift, int bins[]);
int az_check_input__(int data_org[], int options[], double params[],
                     int proc_config[]);
int az_find_index__(int *key, int list[], int *length);
int az_gmax_int__(int *val, int proc_config[]);
int az_gmin_int__(int *val, int proc_config[]);
int az_gsum_int__(int *totals, int proc_config[]);
int az_quick_find__(int *key, int list[],int *length, int *shift, int bins[]);
extern MPI_AZComm *az_get_comm_(int proc_config[]);
extern MPI_AZComm *az_get_comm__(int proc_config[]);


int AZ_using_fortran = AZ_FALSE;

void az_broadcast_(char *ptr, int *length, int proc_config[],int *action)
{
  AZ_broadcast(ptr, *length, proc_config, *action);
}

int  az_check_input_(int data_org[], int options[], double params[],
                     int proc_config[])
{
  /* untouched */
  return(AZ_check_input(data_org, options, params, proc_config));
}

void az_check_msr_(int *bindx, int *N_update, int *N_external,
                   int option, int *proc_config)
{
  AZ_check_msr(bindx, *N_update, *N_external, option, proc_config);
}

void az_check_vbr_(int *N_update, int *N_external, int *option, int bindx[],
                   int bnptr[], int cnptr[], int rnptr[], int proc_config[])
{
  AZ_check_vbr(*N_update, *N_external, *option, bindx, bnptr, cnptr, rnptr,
               proc_config);
}

void az_defaults_(int options[], double params[])
{
  /* untouched */
  AZ_defaults(options, params);
}

void az_exchange_bdry_(double x[], int data_org[], int proc_config[])
{
  /* untouched */
  AZ_exchange_bdry(x, data_org, proc_config);
}

int  az_find_index_(int *key, int list[], int *length)
{
  return(AZ_find_index(*key, list, *length));
}

void az_find_local_indices_(int *N_update, int bindx[], int update[],
                            int *external, int *N_external, int *mat_type,
                            int bpntr[])
{
  AZ_using_fortran = AZ_TRUE;
  AZ_find_local_indices(*N_update, bindx, update, &external,
                        N_external, *mat_type, bpntr);
  AZ_using_fortran = AZ_FALSE;
}

void az_find_procs_for_externs_(int *N_update, int update[], int external[],
                                int *N_external, int proc_config[],
                                int *extern_proc)
{

  AZ_using_fortran = AZ_TRUE;
  AZ_find_procs_for_externs(*N_update, update, external, *N_external,
                            proc_config, &extern_proc);
  AZ_using_fortran = AZ_FALSE;
}

void az_free_memory_(int *label)
{
  AZ_free_memory(*label);
}

double az_gavg_double_(double *var, int proc_config[])
{
  return(AZ_gavg_double(*var, proc_config));
}

double az_gdot_(int *N, double r[], double z[],int proc_config[])
{
  return(AZ_gdot(*N, r, z,proc_config));
}

double az_gmax_double_(double *var, int proc_config[])
{
  return(AZ_gmax_double(*var, proc_config));
}

int az_gmax_int_(int *val, int proc_config[])
{
  return(AZ_gmax_int(*val, proc_config));
}

double az_gmax_matrix_norm_(double val[], int indx[], int bindx[], int rpntr[],
                            int cpntr[], int bpntr[], int proc_config[],
                            int data_org[])
{
  /* unchanged */
  return(AZ_gmax_matrix_norm(val, indx, bindx, rpntr, cpntr, bpntr,
                             proc_config, data_org));
}

double az_gmax_vec_(int *N, double vec[], int proc_config[])
{
  return(AZ_gmax_vec(*N, vec, proc_config));
}

double az_gmin_double_(double *var, int proc_config[])
{
  return(AZ_gmin_double(*var, proc_config));
}

int az_gmin_int_(int *val, int proc_config[])
{
  return(AZ_gmin_int(*val, proc_config));
}

double az_gsum_double_(double *var, int proc_config[])
{
  return(AZ_gsum_double(*var, proc_config));
}

int    az_gsum_int_(int *totals, int proc_config[])
{
  return(AZ_gsum_int(*totals, proc_config));
}

void   az_gsum_vec_int_(int vals[], int vals2[], int *length, int proc_config[])
{
  AZ_gsum_vec_int(vals, vals2, *length, proc_config);
}

double az_gvector_norm_(int *n, int *p, double *x, int *proc_config)
{
  return(AZ_gvector_norm(*n, *p, x, proc_config));
}

void   az_init_quick_find_(int list[], int *length, int *shift, int *bins)
{
  AZ_init_quick_find(list, *length, shift, bins);
}

void az_invorder_vec_(double vector[], int data_org[], int update_index[],
	int rpntr[], double newvec[])
{
  AZ_invorder_vec(vector, data_org, update_index, rpntr, newvec);
}

void   az_matvec_mult_(double *val, int *indx, int *bindx, int *rpntr,
                       int *cpntr, int *bpntr, double *b, double *c,
                       int *exchange_flag, int *data_org)
{
  AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, b, c,
                 *exchange_flag, data_org);
}

void   az_msr2vbr_(double val[], int indx[], int rnptr[], int cnptr[],
                   int bnptr[], int bindx[], int msr_bindx[], double msr_val[],
                   int *total_blk_rows, int *total_blk_cols, int *blk_space,
                   int *nz_space, int *blk_type)
{
  AZ_msr2vbr(val, indx, rnptr, cnptr, bnptr, bindx, msr_bindx, msr_val,
             *total_blk_rows, *total_blk_cols, *blk_space, *nz_space,
             *blk_type);
}

void   az_order_ele_(int update_index[], int extern_index[], int *internal,
                     int *border, int *N_update, int msr_bindx[], int bindx[],
                     int extern_proc[], int *N_external, int *option,
                     int *m_type)
{
  AZ_order_ele(update_index, extern_index, internal, border, *N_update,
               msr_bindx, bindx, extern_proc, *N_external, *option, *m_type);
}

void   az_print_error_(int *error_code)
{
  AZ_print_error(*error_code);
}

void   az_processor_info_(int proc_config[])
{
  /* no change */
  AZ_processor_info(proc_config);
}

int    az_quick_find_(int *key, int list[],int *length, int *shift, int bins[])
{
  return(AZ_quick_find(*key, list,*length, *shift, bins));
}

void   az_read_msr_matrix_(int update[], double *val, int *bindx, int *N_update,
                           int proc_config[])
{
  AZ_using_fortran = AZ_TRUE;
  AZ_read_msr_matrix(update, &val, &bindx, *N_update, proc_config);
  AZ_using_fortran = AZ_FALSE;
}

void   az_read_update_(int *N_update_blks, int update_blks[], int proc_config[],
                       int *bigN, int *chunk, int *input_option)
{
  AZ_using_fortran = AZ_TRUE;
  AZ_read_update(N_update_blks, &update_blks, proc_config, *bigN, *chunk,
                 *input_option);
  AZ_using_fortran = AZ_FALSE;
}

void   az_reorder_matrix_(int *N_update, int bindx[], double val[],
                          int update_index[], int extern_index[], int indx[],
                          int rnptr[], int bnptr[], int *N_external,
                          int cnptr[], int *option, int *mat_type)
{
  AZ_reorder_matrix(*N_update, bindx, val, update_index, extern_index, indx,
		    rnptr, bnptr, *N_external, cnptr, *option, *mat_type);
}

void az_reorder_vec_(double vector[], int data_org[], int update_index[],
        int rpntr[])
{
  AZ_reorder_vec(vector, data_org, update_index, rpntr);
}

MPI_AZComm *az_get_comm_(int proc_config[]) { return(AZ_get_comm(proc_config));}
void az_set_comm_(int proc_config[], MPI_AZComm comm) {
   AZ_set_comm(proc_config, comm); }
void az_set_proc_config_(int proc_config[], MPI_AZComm *comm) {
   AZ_set_proc_config(proc_config, *comm); }

void   az_set_message_info_(int *N_external, int extern_index[],
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

void az_solve_(double x[], double b[], int options[], double params[],
               int indx[], int bindx[], int rpntr[], int cpntr[], int bpntr[],
               double val[], int data_org[], double status[], int proc_config[])
{
  /* no change */
  AZ_solve(x, b, options, params, indx, bindx, rpntr, cpntr, bpntr,
           val, data_org, status, proc_config);
}

void   az_sort_(int list[], int *N, int list2[], double list3[])
{
  AZ_sort(list, *N, list2, list3); 
}

void   az_transform_(int proc_config[],int external[],int bindx[], double val[],
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

void az_broadcast__(char *ptr, int *length, int proc_config[],int *action)
{
  AZ_broadcast(ptr, *length, proc_config, *action);
}

int  az_check_input__(int data_org[], int options[], double params[],
                     int proc_config[])
{
  /* untouched */
  return(AZ_check_input(data_org, options, params, proc_config));
}

void az_check_msr__(int *bindx, int *N_update, int *N_external,
                   int option, int *proc_config)
{
  AZ_check_msr(bindx, *N_update, *N_external, option, proc_config);
}

void az_check_vbr__(int *N_update, int *N_external, int *option, int bindx[],
                   int bnptr[], int cnptr[], int rnptr[], int proc_config[])
{
  AZ_check_vbr(*N_update, *N_external, *option, bindx, bnptr, cnptr, rnptr,
               proc_config);
}

void az_defaults__(int options[], double params[])
{
  /* untouched */
  AZ_defaults(options, params);
}

void az_exchange_bdry__(double x[], int data_org[], int proc_config[])
{
  /* untouched */
  AZ_exchange_bdry(x, data_org, proc_config);
}

int  az_find_index__(int *key, int list[], int *length)
{
  return(AZ_find_index(*key, list, *length));
}

void az_find_local_indices__(int *N_update, int bindx[], int update[],
                            int *external, int *N_external, int *mat_type,
                            int bpntr[])
{
  AZ_using_fortran = AZ_TRUE;
  AZ_find_local_indices(*N_update, bindx, update, &external,
                        N_external, *mat_type, bpntr);
  AZ_using_fortran = AZ_FALSE;
}

void az_find_procs_for_externs__(int *N_update, int update[], int external[],
                                int *N_external, int proc_config[],
                                int *extern_proc)
{

  AZ_using_fortran = AZ_TRUE;
  AZ_find_procs_for_externs(*N_update, update, external, *N_external,
                            proc_config, &extern_proc);
  AZ_using_fortran = AZ_FALSE;
}

void az_free_memory__(int *label)
{
  AZ_free_memory(*label);
}

double az_gavg_double__(double *var, int proc_config[])
{
  return(AZ_gavg_double(*var, proc_config));
}

double az_gdot__(int *N, double r[], double z[],int proc_config[])
{
  return(AZ_gdot(*N, r, z,proc_config));
}

double az_gmax_double__(double *var, int proc_config[])
{
  return(AZ_gmax_double(*var, proc_config));
}

int az_gmax_int__(int *val, int proc_config[])
{
  return(AZ_gmax_int(*val, proc_config));
}

double az_gmax_matrix_norm__(double val[], int indx[], int bindx[], int rpntr[],
                            int cpntr[], int bpntr[], int proc_config[],
                            int data_org[])
{
  /* unchanged */
  return(AZ_gmax_matrix_norm(val, indx, bindx, rpntr, cpntr, bpntr,
                             proc_config, data_org));
}

double az_gmax_vec__(int *N, double vec[], int proc_config[])
{
  return(AZ_gmax_vec(*N, vec, proc_config));
}

double az_gmin_double__(double *var, int proc_config[])
{
  return(AZ_gmin_double(*var, proc_config));
}

int az_gmin_int__(int *val, int proc_config[])
{
  return(AZ_gmin_int(*val, proc_config));
}

double az_gsum_double__(double *var, int proc_config[])
{
  return(AZ_gsum_double(*var, proc_config));
}

int    az_gsum_int__(int *totals, int proc_config[])
{
  return(AZ_gsum_int(*totals, proc_config));
}

void   az_gsum_vec_int__(int vals[], int vals2[],int *length, int proc_config[])
{
  AZ_gsum_vec_int(vals, vals2, *length, proc_config);
}

double az_gvector_norm__(int *n, int *p, double *x, int *proc_config)
{
  return(AZ_gvector_norm(*n, *p, x, proc_config));
}

void   az_init_quick_find__(int list[], int *length, int *shift, int *bins)
{
  AZ_init_quick_find(list, *length, shift, bins);
}

void az_invorder_vec__(double vector[], int data_org[], int update_index[],
	int rpntr[], double newvec[])
{
  AZ_invorder_vec(vector, data_org, update_index, rpntr, newvec);
}

void   az_matvec_mult__(double *val, int *indx, int *bindx, int *rpntr,
                       int *cpntr, int *bpntr, double *b, double *c,
                       int *exchange_flag, int *data_org)
{
  AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, b, c,
                 *exchange_flag, data_org);
}

void   az_msr2vbr__(double val[], int indx[], int rnptr[], int cnptr[],
                   int bnptr[], int bindx[], int msr_bindx[], double msr_val[],
                   int *total_blk_rows, int *total_blk_cols, int *blk_space,
                   int *nz_space, int *blk_type)
{
  AZ_msr2vbr(val, indx, rnptr, cnptr, bnptr, bindx, msr_bindx, msr_val,
             *total_blk_rows, *total_blk_cols, *blk_space, *nz_space,
             *blk_type);
}

void   az_order_ele__(int update_index[], int extern_index[], int *internal,
                     int *border, int *N_update, int msr_bindx[], int bindx[],
                     int extern_proc[], int *N_external, int *option,
                     int *m_type)
{
  AZ_order_ele(update_index, extern_index, internal, border, *N_update,
               msr_bindx, bindx, extern_proc, *N_external, *option, *m_type);
}

void   az_print_error__(int *error_code)
{
  AZ_print_error(*error_code);
}

void   az_processor_info__(int proc_config[])
{
  /* no change */
  AZ_processor_info(proc_config);
}

int    az_quick_find__(int *key, int list[],int *length, int *shift, int bins[])
{
  return(AZ_quick_find(*key, list,*length, *shift, bins));
}

void   az_read_msr_matrix__(int update[],double *val, int *bindx, int *N_update,
                           int proc_config[])
{
  AZ_using_fortran = AZ_TRUE;
  AZ_read_msr_matrix(update, &val, &bindx, *N_update, proc_config);
  AZ_using_fortran = AZ_FALSE;
}

void   az_read_update__(int *N_update_blks,int update_blks[], int proc_config[],
                       int *bigN, int *chunk, int *input_option)
{
  AZ_using_fortran = AZ_TRUE;
  AZ_read_update(N_update_blks, &update_blks, proc_config, *bigN, *chunk,
                 *input_option);
  AZ_using_fortran = AZ_FALSE;
}

void   az_reorder_matrix__(int *N_update, int bindx[], double val[],
                          int update_index[], int extern_index[], int indx[],
                          int rnptr[], int bnptr[], int *N_external,
                          int cnptr[], int *option, int *mat_type)
{
  AZ_reorder_matrix(*N_update, bindx, val, update_index, extern_index, indx,
		    rnptr, bnptr, *N_external, cnptr, *option, *mat_type);
}

void az_reorder_vec__(double vector[], int data_org[], int update_index[],
        int rpntr[])
{
  AZ_reorder_vec(vector, data_org, update_index, rpntr);
}

MPI_AZComm *az_get_comm__(int proc_config[]) {return(AZ_get_comm(proc_config));}
void az_set_comm__(int proc_config[], MPI_AZComm comm) {
   AZ_set_comm(proc_config, comm); }
void az_set_proc_config__(int proc_config[], MPI_AZComm *comm) {
   AZ_set_proc_config(proc_config, *comm); }


void   az_set_message_info__(int *N_external, int extern_index[],
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

void az_solve__(double x[], double b[], int options[], double params[],
               int indx[], int bindx[], int rpntr[], int cpntr[], int bpntr[],
               double val[], int data_org[], double status[], int proc_config[])
{
  /* no change */
  AZ_solve(x, b, options, params, indx, bindx, rpntr, cpntr, bpntr,
           val, data_org, status, proc_config);
}

void   az_sort__(int list[], int *N, int list2[], double list3[])
{
  AZ_sort(list, *N, list2, list3); 
}

void   az_transform__(int proc_config[],int external[],int bindx[], double val[],
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
