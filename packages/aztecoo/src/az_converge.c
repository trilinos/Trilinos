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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "az_aztec.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_compute_global_scalars(AZ_MATRIX *Amat, 
                               double x[], double b[], double r[], double w[],
                               double *r_norm, double *scaled_r_norm,
                               int options[], int data_org[], int proc_config[],
                               int *r_avail, double v1[], double v2[],
                               double *value, 
			       struct AZ_CONVERGE_STRUCT *conv_info)

/*******************************************************************************

  Routine to check against 'eps' for convergence. The type of norm use is
  determined by the variable 'conv_flag' as follows:

                0: ||r||2 / ||r0||2  < eps
                1: ||r||2 / ||b||2   < eps
                2: ||r||2 / ||A||inf < eps
                3: ||r||inf / (||A||inf * ||x||1 + ||b||inf)
                4: ||r/w||2

                where ||*||2 is the Euclidean norm, ||*||inf is the
                infinity norm and ||*|| is the sum norm.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   Aztec User's Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  x:               The current solution vector.

  b:               Right hand side of linear system.

  r:               The current residual vector.

  w:               Weighting vector for convergence norm #4.

  r_norm:          Norm of residual.

  scaled_r_norm:   Norm of residual scaled by norm of the rhs.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  r_avail:         In general, this variable indicates whether or not the
                   residual is available or needs to be made available.
                   In particular,
                      first_time == TRUE     : real residual is available. The
                                               norm of this residual will be
                                               computed and stored in r_norm.
                      first_time == FALSE &&
                      r_avail    == TRUE     : real residual is available. The
                                               norm of this residual will be
                                               computed and stored in r_norm.
                      first_time == FALSE &&
                      r_avail    == FALSE    : real residual is not available.
                                               The norm of the residual is not
                                               computed. Instead, it is assumed
                                               that r_norm is an estimate of the
                                               residual norm.
                   All of this is done for gmres() and tfqmr() where we often
                   have estimates of the residual 2-norm without actually having
                   computed the residual.

                   IMPORTANT: if a convergence test requires a residual norm
                   other than the 2-norm, it is important that
                   AZ_compute_global_scalars() sets r_avail to TRUE. This tells
                   the iterative method (in particular gmres and tfqmr) that the
                   real residual must be computed (at additional cost)
                   and passed in to AZ_compute_global_scalars().

  v1,v2,value:     If v1 != NULL, *value = <v1,v2> where <.,.> is the standard
                   inner product, v1 and v2 are double precision vectors of
                   length data_org[AZ_N_internal] + data_org[AZ_N_border].

                   This is used so that 1 inner product can be grouped together
                   with the inner products required for convergence (to save
                   some messages).

  first_time:      Flag set AZ_TRUE if this is the first time this routine has
                   been called.  Set AZ_FALSE otherwise. See comments for
                   r_avail above for more information.

*******************************************************************************/

{

  /* local variables */

  register int  i;
  static double *temp, *tr;
  static int    total_N;
  int           N, j;
  double        dots[5], tmp[5], dmax, dtemp;
  int           count = 0, one = 1;

  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];
  tr = r;

  if (options[AZ_ignore_scaling]) {
     if ( (conv_info->scaling->action == AZ_left_scaling) || 
          (conv_info->scaling->action == AZ_left_and_right_scaling) ) {
        if (!(*r_avail) && (conv_info->not_initialized==AZ_FALSE)) {
           printf("AZ_compute_global_scalars: Error residual is needed to \
		   ignore scaling in convergence tests\n");
           exit(1);
        }
        *r_avail = AZ_TRUE;
        tr = AZ_manage_memory(N*sizeof(double),AZ_ALLOC,AZ_SYS, "trinconv",&j);
        for (i = 0; i < N; i++) tr[i] = r[i];
        AZ_scale_f(AZ_INVSCALE_RHS, Amat, options, tr, x, proc_config, 
                   conv_info->scaling);
     }
  }

  if (v1 != NULL) dots[count++] = ddot_(&N, v1, &one, v2, &one);

  /* initialize */

  switch (options[AZ_conv]) {

  case AZ_noscaled:
    if ((*r_avail) || conv_info->not_initialized) {
      dots[count] = ddot_(&N, tr, &one, tr, &one);
      AZ_gdot_vec(count + 1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
    }
    else if (v1 != NULL) AZ_gdot_vec(1, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    *scaled_r_norm = *r_norm ;
    break;


  case AZ_r0:
    if ((*r_avail) || conv_info->not_initialized) {
      dots[count] = ddot_(&N, tr, &one, tr, &one);
      AZ_gdot_vec(count + 1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
    }
    else if (v1 != NULL) AZ_gdot_vec(1, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    if (conv_info->not_initialized) {
      conv_info->r0_norm = *r_norm;
      if (fabs(conv_info->r0_norm) < DBL_MIN) conv_info->r0_norm = 1.0;
    }

    *scaled_r_norm = *r_norm / conv_info->r0_norm;
    break;

  case AZ_rhs:
    if (conv_info->not_initialized) {
      dots[count++] = ddot_(&N, tr, &one, tr, &one);

      if  ( (options[AZ_ignore_scaling]) && (tr != r) ) {
          for (i = 0; i < N; i++) tr[i] = b[i];
          AZ_scale_f(AZ_INVSCALE_RHS, Amat, options, tr, x, proc_config, 
                     conv_info->scaling);
      }
      else tr = b;
      dots[count  ] = ddot_(&N, tr, &one, tr, &one);
      AZ_gdot_vec(count + 1, dots, tmp, proc_config);

      conv_info->b_norm = sqrt(dots[count--]);
      if (conv_info->b_norm == 0.0) {
	if ((proc_config[AZ_node]==0) && (options[AZ_output] != AZ_none)) {
	    printf("AZ_compute_global_scalars: ||rhs|| = 0. Can not use AZ_rhs as a convergence option.\n");
            printf("Changing convergence criteria to use unscaled residual norm in convergence tests.\n");
         }
	 conv_info->b_norm = 1.;
      }
      *r_norm = sqrt(dots[count]); 
    }

    else if (*r_avail) {
      dots[count] = ddot_(&N, tr, &one, tr, &one);

      AZ_gdot_vec(count + 1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
    }

    else if (v1 != NULL) AZ_gdot_vec(1, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    *scaled_r_norm = *r_norm / conv_info->b_norm;
    break;

  case AZ_Anorm:
    if (conv_info->not_initialized) {
      dots[count] = ddot_(&N, tr, &one, tr, &one);

      AZ_gdot_vec(count+1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
      if  ((options[AZ_ignore_scaling]) && (conv_info->scaling->A_norm == 0.0)&&
           (options[AZ_scaling] != AZ_none) &&
           (options[AZ_pre_calc] == AZ_reuse)) {
         if ((proc_config[AZ_node]==0) && (options[AZ_output] != AZ_none)) {
            printf("Warning:No previous definition for A_norm found. Was ");
            printf("AZ_iterate used\n\tpreviously and was the scaling object ");
            printf("passed in the same as for\n\tthis invokation of ");
            printf("AZ_iterate()?\n");
         }
       }


      if  ((options[AZ_ignore_scaling]) && (conv_info->scaling->A_norm != 0.0))
         conv_info->A_norm = conv_info->scaling->A_norm;
      else conv_info->A_norm = AZ_gmax_matrix_norm(Amat->val, Amat->indx, 
				   Amat->bindx, Amat->rpntr, Amat->cpntr, 
				   Amat->bpntr, proc_config, Amat->data_org);

      if (fabs(conv_info->A_norm) < DBL_MIN) {
        AZ_p_error("Error: ||A||_infinity = 0\n\tSomething wrong with A\n",
                   proc_config[AZ_node]);
        exit(-1);
      }
    }

    else if (*r_avail) {
      dots[count] = ddot_(&N, tr, &one, tr, &one);

      AZ_gdot_vec(count + 1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
    }

    else if (v1 != NULL) AZ_gdot_vec(1, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    *scaled_r_norm = *r_norm / conv_info->A_norm;
    break;

  case AZ_sol:
    *r_avail = AZ_TRUE;

    if (v1 != NULL) {
      AZ_gdot_vec(1, dots, tmp, proc_config);
      *value = dots[0];
    }

    *r_norm = AZ_gvector_norm(N, -1, tr, proc_config);
    if (conv_info->not_initialized) {
      if  ((options[AZ_ignore_scaling]) && (conv_info->scaling->A_norm == 0.0)&&
           (options[AZ_scaling] != AZ_none) &&
           (options[AZ_pre_calc] == AZ_reuse)) {
         if ((proc_config[AZ_node]==0) && (options[AZ_output] != AZ_none)) {
            printf("Warning:No previous definition for A_norm found. Was ");
            printf("AZ_iterate used\n\tpreviously and was the scaling object ");
            printf("passed in the same as for\n\tthis invokation of ");
            printf("AZ_iterate()?\n");
         }
       }
      if  ((options[AZ_ignore_scaling]) && (conv_info->scaling->A_norm != 0.0))
         conv_info->A_norm = conv_info->scaling->A_norm;
      else conv_info->A_norm = AZ_gmax_matrix_norm(Amat->val, Amat->indx, 
				   Amat->bindx, Amat->rpntr, Amat->cpntr, 
				   Amat->bpntr, proc_config, Amat->data_org);

      if (fabs(conv_info->A_norm) < DBL_MIN) {
        AZ_p_error("Error: ||A||_infinity = 0\n\tSomething wrong with A\n",
                   proc_config[AZ_node]);
        exit(-1);
      }


      if  ( (options[AZ_ignore_scaling]) && (tr != r) ) {
          for (i = 0; i < N; i++) tr[i] = b[i];
          AZ_scale_f(AZ_INVSCALE_RHS, Amat, options, tr, x, proc_config, 
                     conv_info->scaling);
      }
      else tr = b;
      conv_info->b_norm = AZ_gvector_norm(N, -1, tr, proc_config);
    }
    
    if ((options[AZ_ignore_scaling]) && 
        (conv_info->scaling->action == AZ_left_and_right_scaling)) {
          for (i = 0; i < N; i++) tr[i] = x[i];
          AZ_scale_f(AZ_INVSCALE_SOL, Amat, options, b, tr, proc_config, 
                     conv_info->scaling);
    }
    else tr =x;
    *scaled_r_norm = *r_norm / (conv_info->A_norm *
                      AZ_gvector_norm(N, 1, tr, proc_config)+conv_info->b_norm);
    break;

  case AZ_expected_values:
    *r_avail = AZ_TRUE;

    dmax = 0.0; 
    for (i = 0; i < N; i++) {
       dtemp = fabs(tr[i]/w[i]);
       if (dtemp > dmax) dmax = dtemp;
    }
    *scaled_r_norm = AZ_gmax_double(dmax, proc_config);

    if (count != 0) AZ_gdot_vec(count, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    *r_norm        = *scaled_r_norm;
    break;
  case AZ_weighted:
    *r_avail = AZ_TRUE;

    temp = AZ_manage_memory((N + data_org[AZ_N_external]) * sizeof(double),
                             AZ_ALLOC, AZ_SYS,
                             "temp in AZ_compute_global_scalars", &j);
    if (conv_info->not_initialized) {
      total_N = AZ_gsum_int(N, proc_config);
      conv_info->total_N = (double) total_N;
    }
    total_N = conv_info->total_N;

    for (i = 0; i < N; i++) temp[i] = tr[i] / w[i];

    dots[count] = ddot_(&N, temp, &one, temp, &one);

    AZ_gdot_vec(count + 1, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    *scaled_r_norm = sqrt(dots[count] / (double) total_N);
    *r_norm        = *scaled_r_norm;
    break;

  case AZ_inf_noscaled:

    *r_avail = AZ_TRUE;

    if (v1 != NULL) {
      AZ_gdot_vec(1, dots, tmp, proc_config);
      *value = dots[0];
    }

    *r_norm = AZ_gvector_norm(N, -1, r, proc_config);
    *scaled_r_norm = *r_norm;

    break;

  default:
    if (proc_config[AZ_node] == 0) {
      (void) fprintf(stderr, "Error: Improper value, options[AZ_conv] = %d\n",
                     options[AZ_conv]);
    }
    exit(-1);
  }
  conv_info->not_initialized = AZ_FALSE;

} /* AZ_compute_global_scalars */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_scale_true_residual( double x[], double b[],
                            double v[], double w[], double *actual_residual,
                            double *scaled_r_norm, int options[],
                            int data_org[], int proc_config[],
                            AZ_MATRIX *Amat, 
			    struct AZ_CONVERGE_STRUCT *conv_info)


/*******************************************************************************

  Routine to check the true residual and decide whether or not actual
  convergence has been achieved based upon some initial convergence indicators.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   Aztec User's Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  x:               The current solution vector.

  b:               Right hand side of linear system.

  v:

  w:               Weighting vector for convergence norm #4.

  actual_residual: Norm of residual.

  scaled_r_norm:   Norm of residual scaled by norm of the rhs.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int r_avail = AZ_TRUE;

  /**************************** execution begins ******************************/

  /* calculate true residual */

  AZ_compute_residual(b, x, v, proc_config, Amat);

  /* compute scaled residual */

  AZ_compute_global_scalars(Amat, x, b, v, w,
                            actual_residual, scaled_r_norm, options, data_org,
                            proc_config,&r_avail, (double *) 0, (double *) 0,
                            (double *) 0, conv_info);

} /* AZ_scale_true_residual */

struct AZ_CONVERGE_STRUCT *AZ_converge_create()
{
   struct AZ_CONVERGE_STRUCT *temp;

   temp = (struct AZ_CONVERGE_STRUCT *) AZ_allocate(sizeof(struct 
                                                    AZ_CONVERGE_STRUCT));
   if (temp == NULL) {
      printf("AZ_converge_create: Not enough space\n");
      exit(1);
   }
   temp->r0_norm = 0.0;
   temp->A_norm  = 0.0;
   temp->b_norm  = 0.0;
   temp->scaling = NULL;
   temp->total_N       = 0;
   temp->not_initialized = AZ_TRUE;
   return temp;
}
void AZ_converge_destroy(struct AZ_CONVERGE_STRUCT **temp)
{
   if (*temp != NULL) AZ_free(*temp);
   *temp = NULL;
}
