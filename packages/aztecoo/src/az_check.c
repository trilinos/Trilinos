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
#ifndef lint
static char rcsid[] = "$Id$";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "az_aztec.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_check_input(int data_org[], int options[], double params[],
                   int proc_config[])

/*******************************************************************************

  Routine to perform checks for iterative solver library. This is to be called
  by the user of the solver library who must supply the necessary information
  in the input arrays. The routine checks that these values. If all the values
  are valid AZ_check_input() returns 0. Otherwise it returns an error code.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     int, error code, 0 => no errors.
  ============

  Parameter list:
  ===============

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int      i, sum;
  char     yo[32];

  sprintf(yo, "AZ_check_input: ");

  /* set a few default values */

  if (params[AZ_tol]            < 0.0       ) params[AZ_tol]          = 1.e-06;
  if (params[AZ_drop]           < 0.0       ) params[AZ_drop]         = 0.;
  if (params[AZ_omega]< 0.0 || params[AZ_omega]>1.) 
                                              params[AZ_omega]        = 1.;
  if (data_org[AZ_N_border]    == AZ_default) data_org[AZ_N_border]   = 0;
  if (data_org[AZ_N_external]  == AZ_default) data_org[AZ_N_external] = 0;
  if (data_org[AZ_N_bord_blk]  == AZ_default) data_org[AZ_N_bord_blk] = 0;
  if (data_org[AZ_N_ext_blk]   == AZ_default) data_org[AZ_N_ext_blk]  = 0;
  if (data_org[AZ_N_neigh]     == AZ_default) data_org[AZ_N_neigh]    = 0;
  if (data_org[AZ_total_send]  == AZ_default) data_org[AZ_total_send] = 0;
  if (data_org[AZ_name]        == AZ_default) data_org[AZ_name]       = 1;
  if (data_org[AZ_matrix_type] == AZ_default) data_org[AZ_matrix_type] =
                                                AZ_VBR_MATRIX;

  sum = 0;
  for (i = 0; i < data_org[AZ_N_neigh]; i++)
    sum += data_org[AZ_send_length + i];
  data_org[AZ_total_send] = sum;

  /* check for warnings */

  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
    if ((data_org[AZ_N_int_blk]  != data_org[AZ_N_internal])  ||
        (data_org[AZ_N_bord_blk] != data_org[AZ_N_border]) ||
        (data_org[AZ_N_ext_blk]  != data_org[AZ_N_external])) {

      /* set blocks information for msr applications */

      (void) fprintf(stdout, "Warning: Setting the number of blocks equal\n");
      (void) fprintf(stdout, "         to the number of unknowns for this\n");
      (void) fprintf(stdout, "         MSR application.                  \n");
      data_org[AZ_N_int_blk]  = data_org[AZ_N_internal];
      data_org[AZ_N_bord_blk] = data_org[AZ_N_border];
      data_org[AZ_N_ext_blk]  = data_org[AZ_N_external];
    }
  }

  /* do some error checking on the contents of "options" */

  if ( options[AZ_solver]     < AZ_cg           ||
       options[AZ_solver]     > AZ_lu            ) return  -1;

  if (options[AZ_scaling]     < AZ_none         ||
      options[AZ_scaling]     > AZ_sym_BJacobi   ) return  -2;

  if ((options[AZ_precond]    < AZ_none         ||
      options[AZ_precond]     > AZ_user_precond  ) &&
     (options[AZ_precond]    >= AZ_SOLVER_PARAMS)) return  -3;

  if (options[AZ_conv]        < AZ_r0                 ||
      options[AZ_conv]        > AZ_inf_noscaled  ) return  -4;

  if (options[AZ_output]      < AZ_all           ) return  -5;

  if (options[AZ_pre_calc]    < AZ_recalc       &&
      options[AZ_pre_calc]    > AZ_sys_reuse     ) return  -6;

  if (options[AZ_max_iter]    < 1                ) return  -7;

  if (options[AZ_precond]    == AZ_ls           &&
      options[AZ_poly_ord]    > AZ_MAX_POLY_ORDER) return  -8;

  if (options[AZ_overlap]  < AZ_diag )             return -9;

  if (options[AZ_solver]     == AZ_gmres        &&
      options[AZ_kspace]      < 1                ) return -10;

  if (options[AZ_solver]     == AZ_GMRESR       &&
      options[AZ_kspace]      < 1                ) return -10;

  if (options[AZ_orthog]     != AZ_classic      &&
      options[AZ_orthog]     != AZ_modified      ) return -11;

  if (options[AZ_aux_vec]    != AZ_resid        &&
      options[AZ_aux_vec]    != AZ_rand          ) return -12;

  if (data_org[AZ_N_border]   < 0                ) return -13;
  if (data_org[AZ_N_internal] < 0                ) return -14;
  if (data_org[AZ_N_external] < 0                ) return -15;
  if (data_org[AZ_N_bord_blk] < 0                ) return -16;
  if (data_org[AZ_N_int_blk]  < 0                ) return -17;
  if (data_org[AZ_N_ext_blk]  < 0                ) return -18;

  if (data_org[AZ_N_neigh]    < 0               ||
      data_org[AZ_N_neigh]    > AZ_MAX_NEIGHBORS ) return -19;

  if (proc_config[AZ_N_procs]<= 0                ) return -20;
  if (proc_config[AZ_node]    < 0                ) return -21;

  /* vector of neighboring processor numbers */

  for (i = 0; i < data_org[AZ_N_neigh]; i++) {
    if (data_org[AZ_neighbors + i] >= proc_config[AZ_N_procs] ||
        data_org[AZ_neighbors + i] < 0           ) return -22;
  }

  /* vector of number of unknowns to receive from neighbor */

  for (i = 0; i < data_org[AZ_N_neigh]; i++) {
    if (data_org[AZ_rec_length + i ] > data_org[AZ_N_external] ||
        data_org[AZ_rec_length + i ] < 0         ) return -23;
  }

  /* vector of number of unknowns to send to neighbor */

  sum = data_org[AZ_N_internal] + data_org[AZ_N_border];

  for (i = 0; i < data_org[AZ_N_neigh]; i++) {
    if (data_org[AZ_send_length+i] > sum ||
        data_org[AZ_send_length + i] < 0         ) return -24;
    else {
      if (data_org[AZ_send_length + i] > data_org[AZ_N_border]) {
        (void) fprintf(stderr, "WARNING: Processor %d sends more than just its \
border points implying that the\n         matrix sparsity pattern is not \
symmetric.\n", 
                       proc_config[AZ_node]);
/*
        (void) fprintf(stderr, "WARNING: Processor %d sends more than just ",
                       proc_config[AZ_node]);
        (void) fprintf(stderr,"its border points implying that the matrix\n");
        (void) fprintf(stderr,"         sparsity pattern is not symmetric.\n");
*/
      }
    }
  }

  if ( (options[AZ_output] > 0) && proc_config[AZ_node] == 0) {
    (void) fprintf(stdout,"\n======================================="
                   "========================================\n");
    (void) fprintf(stdout,"%sSetup information on processor 0\n\n", yo);
    (void) fprintf(stdout,"\tsolver:\t\t\t\t\t%d\n", options[AZ_solver]);
    (void) fprintf(stdout,"\tconvergence flag:\t\t\t%d\n", options[AZ_conv]);
    (void) fprintf(stdout,"\tmaximum iterations:\t\t\t%d\n",
                   options[AZ_max_iter]);
    (void) fprintf(stdout,"\treordering:    \t\t\t\t%d\n", options[AZ_reorder]);
    (void) fprintf(stdout,"\tpreconditioner:\t\t\t\t%d\n", options[AZ_precond]);
    (void) fprintf(stdout,"\tpolynomial order:\t\t\t%d\n",
                   options[AZ_poly_ord]);
  if (options[AZ_solver]==AZ_gmres) {
    (void) fprintf(stdout,"\tGMRES ill conditioning threshold:\t\t\t%7.1e\n",
                   params[AZ_ill_cond_thresh]);
    (void) fprintf(stdout,"\tGMRES restart size:\t\t\t%d\n",
                   options[AZ_kspace]);
    (void) fprintf(stdout,"\torthogonalization:\t\t\t%d\n", options[AZ_orthog]);}
    (void) fprintf(stdout,"\ttolerance:\t\t\t\t%7.1e\n", params[AZ_tol]);
    (void) fprintf(stdout,"\tdrop:\t\t\t\t\t%7.1e\n", params[AZ_drop]);
    if ( (options[AZ_precond] == AZ_dom_decomp) &&
         (options[AZ_subdomain_solve] == AZ_ilut) ) {
      if (params[AZ_rthresh]!=0.0)
      (void) fprintf(stdout,"\tRelative threshold:\t\t\t%7.1e\n", params[AZ_rthresh]);
      if (params[AZ_athresh]!=0.0)
      (void) fprintf(stdout,"\tAbsolute threshold:\t\t\t%7.1e\n", params[AZ_athresh]);

    (void) fprintf(stdout,"\tfill-in:\t\t\t\t%7.1e\n", params[AZ_ilut_fill]);}
    if ( (options[AZ_precond] == AZ_dom_decomp) &&
         (options[AZ_subdomain_solve] == AZ_rilu) ) {
    (void) fprintf(stdout,"\tomega:\t\t\t\t\t%7.1e\n", params[AZ_omega]);}
    if ( (options[AZ_precond] == AZ_dom_decomp) && (
         (options[AZ_subdomain_solve] == AZ_rilu) ||
         (options[AZ_subdomain_solve] == AZ_ilu ) ||
         (options[AZ_subdomain_solve] == AZ_bilu) ||
         (options[AZ_subdomain_solve] == AZ_bilu_ifp) ||
         (options[AZ_subdomain_solve] == AZ_icc)) ) {
      if (params[AZ_rthresh]!=0.0)
      (void) fprintf(stdout,"\tRelative threshold:\t\t\t%7.1e\n", params[AZ_rthresh]);
      if (params[AZ_athresh]!=0.0)
      (void) fprintf(stdout,"\tAbsolute threshold:\t\t\t%7.1e\n", params[AZ_athresh]);
    (void) fprintf(stdout,"\tfill-in:\t\t\t\t\t%d\n", options[AZ_graph_fill]);
    (void) fprintf(stdout,"\toverlap:\t\t\t\t\t%d\n", options[AZ_overlap]);
  if ( proc_config[AZ_N_procs] > 0) 
    (void) fprintf(stdout,"\ttypeoverlap:\t\t\t\t%d\n", options[AZ_type_overlap]);}
    (void) fprintf(stdout, "\n");
  }

  /* output debug information */

 if ((options[AZ_output] > 0) && proc_config[AZ_node] == 0) {
    (void) fprintf(stdout, "\tNumber of internal unknowns:\t\t%d\n",
                   data_org[AZ_N_internal]);
    (void) fprintf(stdout, "\tNumber of border  unknowns:\t\t%d\n",
                   data_org[AZ_N_border]);
    (void) fprintf(stdout, "\tTotal number of unknowns:\t\t%d\n",
                   data_org[AZ_N_internal] + data_org[AZ_N_border]);
    (void) fprintf(stdout,"\tNumber of external unknowns:\t\t%d\n",
                   data_org[AZ_N_external]);
    (void) fprintf(stdout, "\tNumber of internal blocks:\t\t%d\n",
                   data_org[AZ_N_int_blk]);
    (void) fprintf(stdout, "\tNumber of border  blocks:\t\t%d\n",
                   data_org[AZ_N_bord_blk]);
    (void) fprintf(stdout, "\tTotal number of blocks:\t\t\t%d\n",
                   data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk]);
    (void) fprintf(stdout,"\tNumber of external blocks:\t\t%d\n",
                   data_org[AZ_N_ext_blk]);
    (void) fprintf(stdout, "\tNumber of processors:\t\t\t%d\n",
                   proc_config[AZ_N_procs]);
    (void) fprintf(stdout, "\tNode number:\t\t\t\t%d\n", proc_config[AZ_node]);
    (void) fprintf(stdout, "\tNumber of neighbors:\t\t\t%d\n",
                   data_org[AZ_N_neigh]);
    (void) fprintf(stdout, "\tNumber of unknowns sent to neighbors:\t%d\n",
                   data_org[AZ_total_send]);
    (void) fprintf(stdout, "======================================="
                   "========================================\n");
  }

  /* no errors detected */

  return 0;

} /* AZ_check_input() */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_print_error(int error_code)

/*******************************************************************************

  Print error message corresponding to 'error_code'. Typically, 'error_code'
  is generated by AZ_check_input().

  Author:          Ray S. Tuminaro, 1422, SNL
  =======

  Parameter list:
  ===============
  error_code:      Value generated by AZ_check_input() which indicates
                   an error in one of the arrays:
                       options, params, data_org, proc_config

*******************************************************************************/

{

  char yo[32];
  char options_str[32];
  char data_org_str[32];
  char proc_config_str[32];
  char value_str[32];

  sprintf(yo, "AZ_print_error: ");
  sprintf(options_str, "invalid options[");
  sprintf(data_org_str, "invalid data_org[");
  sprintf(proc_config_str, "invalid proc_config[");
  sprintf(value_str, "] value");

  switch(error_code) {
  case 0:
    break;
  case -1:
    (void) fprintf(stderr, "%s%sAZ_solver%s\n", yo, options_str, value_str);
  break;
  case -2:
    (void) fprintf(stderr, "%s%sAZ_scaling%s\n", yo, options_str, value_str);
  break;
  case -3:
    (void) fprintf(stderr, "%s%sAZ_precond%s\n", yo, options_str, value_str);
  break;
  case -4:
    (void) fprintf(stderr, "%s%sAZ_conv%s\n", yo, options_str, value_str);
  break;
  case -5:
    (void) fprintf(stderr, "%s%sAZ_output%s\n", yo, options_str, value_str);
  break;
  case -6:
    (void) fprintf(stderr, "%s%sAZ_pre_calc%s\n", yo, options_str, value_str);
  break;
  case -7:
    (void) fprintf(stderr, "%s%sAZ_max_iter%s\n", yo, options_str, value_str);
  break;
  case -8:
    (void) fprintf(stderr, "%s%sAZ_poly_ord%s\n", yo, options_str, value_str);
  break;
  case -9:
    (void) fprintf(stderr, "%s%sAZ_overlap%s\n", yo, options_str, value_str);
  break;
  case -10:
    (void) fprintf(stderr, "%s%sAZ_kspace%s\n", yo, options_str, value_str);
  break;
  case -11:
    (void) fprintf(stderr, "%s%sAZ_orthog%s\n", yo, options_str, value_str);
  break;
  case -12:
    (void) fprintf(stderr, "%s%sAZ_aux_vec%s\n", yo, options_str, value_str);
  break;
  case -13:
    (void) fprintf(stderr, "%s%sAZ_N_border%s\n", yo, data_org_str, value_str);
  break;
  case -14:
    (void) fprintf(stderr, "%s%sAZ_N_internal%s\n", yo, data_org_str,
                   value_str);
  break;
  case -15:
    (void) fprintf(stderr, "%s%sAZ_N_external%s\n", yo, data_org_str,
                   value_str);
  break;
  case -16:
    (void) fprintf(stderr, "%s%sAZ_N_bord_blk%s\n", yo, data_org_str,
                   value_str);
  break;
  case -17:
    (void) fprintf(stderr, "%s%sAZ_N_int_blk%s\n", yo, data_org_str, value_str);
  break;
  case -18:
    (void) fprintf(stderr, "%s%sAZ_N_ext_blk%s\n", yo, data_org_str, value_str);
  break;
  case -19:
    (void) fprintf(stderr, "%s%sAZ_N_neigh%s\n", yo, data_org_str, value_str);
  break;
  case -20:
    (void) fprintf(stderr, "%s%sAZ_N_procs%s\n", yo, proc_config_str,
                   value_str);
  break;
  case -21:
    (void) fprintf(stderr, "%s%sAZ_N_node%s\n", yo, proc_config_str, value_str);
  break;
  case -22:
    (void) fprintf(stderr, "%s%sAZ_neighbors+i%s\n", yo, data_org_str,
                   value_str);
  break;
  case -23:
    (void) fprintf(stderr, "%s%sAZ_rec_length+i%s\n", yo, data_org_str,
                   value_str);
  break;
  case -24:
    (void) fprintf(stderr, "%s%sAZ_send_length+i%s\n", yo, data_org_str,
                   value_str);
  break;
  default:
    (void) fprintf(stderr, "%sERROR: invalid error code (%d)\n", yo,
                   error_code);
  }

} /* AZ_print_error */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_check_options(int options[], int az_proc, int data_org[], int az_nprocs,
                     double params[], AZ_MATRIX *Amat, AZ_PRECOND *precond)

/*******************************************************************************

  Check values in option arrays for validity and compatibility.

  Author:          Ray S. Tuminaro, 1422, SNL
  =======

  Return code:     int, 1 indicates satisfactory check, 0 indicates a warning or
  ============     error has been encountered.

  Parameter list:
  ===============

  options:         Determines specific solution method and other parameters.

  az_proc:         Current processor.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

  az_nprocs:       Number of processor in the current machine configuration.

  params:          Drop tolerance and convergence tolerance info.

*******************************************************************************/

{

  int   i;
  char     yo[32];

  int *sub_options;
  double     *sub_params, *sub_status;
  AZ_MATRIX  *sub_matrix;
  AZ_PRECOND *sub_precond;
  struct AZ_SCALING *sub_scaling;

  /**************************** execution begins ******************************/

  sprintf(yo, "AZ_check_options: ");
  switch (options[AZ_solver]) {

  case AZ_cg:
    break;
  case AZ_analyze:
    break;
  case AZ_GMRESR:
    break;
  case AZ_gmres:
    break;
  case AZ_cgs:
    break;
  case AZ_tfqmr:
    break;
  case AZ_bicgstab:
    break;
  case AZ_symmlq:
    break;
  case AZ_fixed_pt:
    break;
  case AZ_lu:
    if (az_nprocs != 1) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sERROR: LU not implemented in parallel."
                       "\n       Try domain decompostion with LU "
                       "preconditioning.\n\n", yo);
      }
      return 0;
    }

/*
    if ((data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
        (data_org[AZ_matrix_type] != AZ_VBR_MATRIX)) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "ERROR: LU factorization can only be "
                       "used with MSR or VBR matrices.\n"
                       "       data_org[AZ_matrix_type] = %d\n\n",
                       data_org[AZ_matrix_type]);
      }
      return 0;
    }
*/
    break;

  default:
    if (az_proc == 0) {
      (void) fprintf(stderr, "%sERROR: options[AZ_solver] has improper value "
                     "(%d)\n\n", yo, options[AZ_solver]);
    }
    return 0;
  }

  switch (options[AZ_precond]) {

  case AZ_none:
    break;
  case AZ_Neumann:
    break;
  case AZ_ls:
    if (options[AZ_poly_ord] > AZ_MAX_POLY_ORDER) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sERROR: Exceeds highest order least_squares "
                       "polynomial available: %d vs %d\n\n", yo,
                       options[AZ_poly_ord], AZ_MAX_POLY_ORDER);
      }
      return 0;
    }
    break;

  case AZ_Jacobi:
    break;

  case AZ_sym_GS:
    if (data_org[AZ_matrix_type] != AZ_MSR_MATRIX) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sERROR: sym GS preconditioning can only "
                       "be used with MSR matrices.\n"
                       "       data_org[AZ_matrix_type] = %d\n\n", yo,
                       data_org[AZ_matrix_type]);
      }
      return 0;
    }
/*
    if (precond->Pmat!=Amat) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sERROR: sym GS preconditioning can only %s",
                       yo,"be used with Pmat=Amat .\n");
      }
      return 0;
    }
*/



    break;

  case AZ_bilu:
/* Begin Aztec 2.1 mheroux mod */
  case AZ_bilu_ifp:
/* End Aztec 2.1 mheroux mod */
     if (options[AZ_reorder]) {
        options[AZ_reorder] = 0;
        if ((options[AZ_output] != AZ_none) && (az_proc  == 0)) {
           printf("\t\t***** Reordering not implemented for Block ILU.\n");
           printf("\t\t***** Continuing without reordering\n");
        }
     }

    if (data_org[AZ_matrix_type] != AZ_VBR_MATRIX) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sERROR: Block ILU can only be used on VBR "
                       "matrices\n       data_org[AZ_matrix_type] = %d\n", yo,
                       data_org[AZ_matrix_type]);
        (void) fprintf(stderr, "       options[AZ_precond]      = %d\n\n",
                       options[AZ_precond]);
      }
      return 0;
    }
    if ( (options[AZ_solver] == AZ_cg) && (options[AZ_overlap] > 0) &&
         (options[AZ_type_overlap] != AZ_symmetric)) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: Preconditioned matrix may not be"
                       " symmetric (due to overlap).\n\n", yo);
      }
    }
    break;
  case AZ_multilevel:
    if (az_proc == 0) 
      printf("Are you sure you want the multilevel preconditioner\n");
    break;
  case AZ_dom_decomp:
/* Begin Aztec 2.1 mheroux mod */
    if ((options[AZ_subdomain_solve]==AZ_bilu ||
        options[AZ_subdomain_solve]==AZ_bilu_ifp)&&(options[AZ_reorder])){
/* End Aztec 2.1 mheroux mod */
        options[AZ_reorder] = 0;
        if ((options[AZ_output] != AZ_none) && (az_proc  == 0)) {
           printf("\t\t***** Reordering not implemented for Block ILU.\n");
           printf("\t\t***** Continuing without reordering\n");
        }
    }
   if ( (options[AZ_solver] == AZ_cg) && 
          ((options[AZ_subdomain_solve] == AZ_ilu) ||
           (options[AZ_subdomain_solve] == AZ_rilu) ||
           (options[AZ_subdomain_solve] == AZ_bilu_ifp) ||
           (options[AZ_subdomain_solve] == AZ_ilut))) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: Preconditioned matrix may not be"
                       " symmetric.\n\n", yo);
      }
    }

    if ( (options[AZ_subdomain_solve] != AZ_lu  ) &&
         (options[AZ_subdomain_solve] != AZ_ilu ) &&
         (options[AZ_subdomain_solve] != AZ_icc ) &&
         (options[AZ_subdomain_solve] != AZ_rilu) &&
         (options[AZ_subdomain_solve] != AZ_ilut) &&
/* Begin Aztec 2.1 mheroux mod */
         (options[AZ_subdomain_solve] != AZ_bilu_ifp) &&
/* End Aztec 2.1 mheroux mod */
         (options[AZ_subdomain_solve] != AZ_bilu) ) {
       if (options[AZ_subdomain_solve] >= AZ_SOLVER_PARAMS) {
          if (az_proc == 0) {
             (void) fprintf(stderr, "%sERROR: options[AZ_subdomain_solve]"
                    " has improper value = %d\n\n", yo, 
                    options[AZ_subdomain_solve]);
           }
           return 0; 
        }
        else {
           AZ_recover_sol_params(options[AZ_subdomain_solve], 
			         &sub_options,&sub_params,&sub_status,
                                 &sub_matrix,&sub_precond,&sub_scaling);
           if (!AZ_check_options(sub_options, az_proc, data_org, az_nprocs,
				    sub_params, sub_matrix, sub_precond)) 
	      return 0;
        }
    }
#ifdef eigen
  case AZ_slu:
#endif
  case AZ_lu:
  case AZ_ilu:
  case AZ_icc:
  case AZ_rilu:
  case AZ_ilut:
    if (options[AZ_overlap] < AZ_diag) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sERROR: Negative overlapping not allowed\n",
                       yo);
      }
       return 0;
    }
    if ( (options[AZ_solver] == AZ_cg) && (options[AZ_overlap] > 0) &&
         (options[AZ_type_overlap] != AZ_symmetric)) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: Preconditioned matrix may not be"
                       " symmetric (due to overlap).\n\n", yo);
      }
    }

    break;

  case AZ_smoother:
    break;
  case AZ_user_precond:
    break;

  default:
    if (options[AZ_precond] >= AZ_SOLVER_PARAMS) {
    if (az_proc == 0) {
      (void) fprintf(stderr, "%sERROR: options[AZ_precond] has improper "
                     "value = %d\n\n", yo, options[AZ_precond]);
    }
    return 0; }
    else {
      AZ_recover_sol_params(options[AZ_precond], &sub_options, &sub_params,
                            &sub_status, &sub_matrix, &sub_precond, &sub_scaling);
      if (!AZ_check_options(sub_options, az_proc, data_org, az_nprocs,
                     sub_params, sub_matrix, sub_precond)) return 0;
    }
  }

  switch (options[AZ_scaling]) {

  case AZ_none:
    break;
  case AZ_Jacobi:
    if ((data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
        (data_org[AZ_matrix_type] != AZ_VBR_MATRIX)) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: Jacobi scaling can only be "
                       "with MSR or VBR\n                 matrices. "
                       "Turning off scaling.\n\n", yo);
      }
    }

    if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: Jacobi scaling for VBR matrices "
                       "is not implemented. Substituting\n"
                       "         block Jacobi instead.\n\n", yo);
      }
    }

    if (options[AZ_solver] == AZ_cg) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: Jacobi scaling may make matrix "
                       "unsymmetric for CG.\n\n", yo);
      }
    }
    break;

  case AZ_BJacobi:
    if ((data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
        (data_org[AZ_matrix_type] != AZ_VBR_MATRIX)) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: block Jacobi scaling can only be "
                       "used with MSR or VBR\n                 matrices."
                       "Turning off scaling.\n\n", yo);
      }
    }

    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: Block Jacobi for MSR matrices is "
                       "not implemented. Substituting \n"
                       "         Jacobi instead.\n\n", yo);
      }
    }

    if (options[AZ_solver] == AZ_cg) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: Jacobi scaling may make matrix "
                       "unsymmetric for CG.\n\n", yo);
      }
    }
    break;

  case AZ_row_sum:
    if ((data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
        (data_org[AZ_matrix_type] != AZ_VBR_MATRIX)) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: left row-sum scaling can only be "
                       "with MSR or VBR\n                  matrices. "
                       "Turning off.\n\n", yo);
      }
    }

    if (options[AZ_solver] == AZ_cg) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: row sum scaling may make matrix "
                       "unsymmetric for CG.\n\n", yo);
      }
    }
    break;

  case AZ_sym_diag:
    if ((data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
        (data_org[AZ_matrix_type] != AZ_VBR_MATRIX)) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: sym diag scaling can only be "
                       "used with MSR or VBR\n                 matrices. "
                       "Turning off.\n\n", yo);
      }
    }
    break;


  case AZ_sym_BJacobi:
    if (data_org[AZ_matrix_type] != AZ_VBR_MATRIX) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: sym block diag. scaling can only be "
                       "used with VBR\n                 matrices."
                       " Turning off.\n\n", yo);
      }
    }
    break;

  case AZ_sym_row_sum:
    if ((data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
        (data_org[AZ_matrix_type] != AZ_VBR_MATRIX)) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: sym row scaling can only be "
                       "used with MSR or VBR matrices.\n"
                       "                     Turning off.\n\n", yo);
      }
    }
    break;

  case AZ_equil:
    if ((data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
        (data_org[AZ_matrix_type] != AZ_VBR_MATRIX)) {
      if (az_proc == 0) {
        (void) fprintf(stderr, "%sWARNING: equilibrated scaling can only be "
                       "used with MSR or VBR matrices.\n"
                       "                     Turning off.\n\n", yo);
      }
    }
    break;

  default:
    if (az_proc == 0) {
      (void) fprintf(stderr, "%sERROR: scaling flag %d not implemented.\n\n",
                     yo, options[AZ_scaling]);
    }
    return 0;
  }

  /* check the the norm used */
  if (((options[AZ_conv] == AZ_sol) || (options[AZ_conv] == AZ_Anorm)) &&
      ((data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
       (data_org[AZ_matrix_type] != AZ_VBR_MATRIX))) {
        if (az_proc == 0) {
        (void) fprintf(stderr, "%sERROR: The norm of the matrix can only be "
                       "computed  with MSR or VBR matrices.\n"
                       "       data_org[AZ_matrix_type] = %d\n\n", yo,
                       data_org[AZ_matrix_type]);
      }
      return 0;
  }



  if (((options[AZ_conv] == 4) || (options[AZ_conv] == 3)) &&
      ((options[AZ_solver] == AZ_gmres) || (options[AZ_solver] == AZ_tfqmr))){
    if (az_proc == 0) {
      (void) fprintf(stderr, "%sWARNING: This convergence option requires "
                     "slightly more work for this solver\n"
                     "         in order to compute the residual.\n\n", yo);
    }
  }

  /* check the weights */

  if (options[AZ_conv] == 4) {
    for (i = 0; i < data_org[AZ_N_internal] + data_org[AZ_N_border]; i++) {
      if (params[AZ_weights+i] <= 0.0) {
        (void) fprintf(stderr, "%sWARNING: A weight vector component is <= 0, "
                       "check params[AZ_WEIGHTS]\n\n", yo);
        return 0;
      }
    }
  }

  return 1;

} /* AZ_check_options */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_defaults(int options[], double params[])

/*******************************************************************************

  Routine to set default parameters for iterative solver options.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

*******************************************************************************/

{
  int i;

  for (i = 0 ; i < AZ_OPTIONS_SIZE; i++ ) options[i] = 0;
  for (i = 0 ; i < AZ_PARAMS_SIZE; i++ )  params[i] = 0.0;

  /* setup default options */

  options[AZ_solver]   = AZ_gmres;
  options[AZ_scaling]  = AZ_none;
  options[AZ_precond]  = AZ_none;
  options[AZ_conv]     = AZ_r0;
  options[AZ_output]   = 1;
  options[AZ_pre_calc] = AZ_calc;
  options[AZ_max_iter] = 500;
  options[AZ_poly_ord] = 3;
  options[AZ_overlap]  = AZ_none;
  options[AZ_type_overlap]  = AZ_standard;
  options[AZ_kspace]   = 30;
  options[AZ_orthog]   = AZ_classic;
  options[AZ_aux_vec]  = AZ_resid;
  options[AZ_reorder]  = 1;
  options[AZ_keep_info]= 0;
  options[AZ_subdomain_solve] = AZ_ilut;
  options[AZ_graph_fill] = 0;
  options[AZ_init_guess] = AZ_NOT_ZERO;
  options[AZ_keep_kvecs] = 0;
  options[AZ_apply_kvecs]= AZ_FALSE;
  options[AZ_orth_kvecs] = AZ_FALSE;
  options[AZ_ignore_scaling] = AZ_FALSE;
  options[AZ_check_update_size] = AZ_FALSE;
  options[AZ_extreme] = AZ_high;

  params[AZ_tol]  = 1.0e-06;
  params[AZ_drop] = 0.0;
  params[AZ_ilut_fill] = 1.;
  params[AZ_omega]= 1.;
  params[AZ_rthresh]= 0.;
  params[AZ_athresh]= 0.;
  params[AZ_update_reduction] = 10e10;
  params[AZ_ill_cond_thresh] = 1.0e11;
}
