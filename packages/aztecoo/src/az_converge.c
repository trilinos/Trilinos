/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "az_aztec.h"
#include "az_blas_wrappers.h"

extern int az_iterate_id;

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
  int print_info;

  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];
  tr = r;

  if (options[AZ_ignore_scaling]) {
    if ( (conv_info->scaling->action == AZ_left_scaling) ||
         (conv_info->scaling->action == AZ_left_and_right_scaling) ) {
      if (!(*r_avail) && (conv_info->not_initialized==AZ_FALSE)) {
        AZ_printf_err("AZ_compute_global_scalars: Error residual is needed to \
                   ignore scaling in convergence tests\n");
        exit(1);
      }
      *r_avail = AZ_TRUE;
      tr = AZ_manage_memory(N*sizeof(double),AZ_ALLOC,AZ_SYS+az_iterate_id, "trinconv",&j);
      for (i = 0; i < N; i++) tr[i] = r[i];
      AZ_scale_f(AZ_INVSCALE_RHS, Amat, options, tr, x, proc_config,
                 conv_info->scaling);
    }
  }

  if (v1 != NULL) dots[count++] = DDOT_F77(&N, v1, &one, v2, &one);

  /* initialize */

  switch (options[AZ_conv]) {

  case AZ_noscaled:
    if ((*r_avail) || conv_info->not_initialized) {
      dots[count] = DDOT_F77(&N, tr, &one, tr, &one);
      AZ_gdot_vec(count + 1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
    }
    else if (v1 != NULL) AZ_gdot_vec(1, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    *scaled_r_norm = *r_norm ;
    break;


  case AZ_r0:
    if ((*r_avail) || conv_info->not_initialized) {
      dots[count] = DDOT_F77(&N, tr, &one, tr, &one);
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
      dots[count++] = DDOT_F77(&N, tr, &one, tr, &one);

      if  ( (options[AZ_ignore_scaling]) && (tr != r) ) {
        for (i = 0; i < N; i++) tr[i] = b[i];
        AZ_scale_f(AZ_INVSCALE_RHS, Amat, options, tr, x, proc_config,
                   conv_info->scaling);
      }
      else tr = b;
      dots[count  ] = DDOT_F77(&N, tr, &one, tr, &one);
      AZ_gdot_vec(count + 1, dots, tmp, proc_config);

      conv_info->b_norm = sqrt(dots[count--]);
      if (conv_info->b_norm == 0.0) {
        if ((proc_config[AZ_node]==0) && (options[AZ_output] != AZ_none)) {
          AZ_printf_out("AZ_compute_global_scalars: ||rhs|| = 0. Can not use AZ_rhs as a convergence option.\n");
          AZ_printf_out("Changing convergence criteria to use unscaled residual norm in convergence tests.\n");
        }
        conv_info->b_norm = 1.;
      }
      *r_norm = sqrt(dots[count]);
    }

    else if (*r_avail) {
      dots[count] = DDOT_F77(&N, tr, &one, tr, &one);

      AZ_gdot_vec(count + 1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
    }

    else if (v1 != NULL) AZ_gdot_vec(1, dots, tmp, proc_config);

    if (v1 != NULL) *value = dots[0];

    *scaled_r_norm = *r_norm / conv_info->b_norm;
    break;

  case AZ_Anorm:
    if (conv_info->not_initialized) {
      dots[count] = DDOT_F77(&N, tr, &one, tr, &one);

      AZ_gdot_vec(count+1, dots, tmp, proc_config);
      *r_norm = sqrt(dots[count]);
      if  ((options[AZ_ignore_scaling]) && (conv_info->scaling->A_norm == 0.0)&&
           (options[AZ_scaling] != AZ_none) &&
           (options[AZ_pre_calc] == AZ_reuse)) {
        if ((proc_config[AZ_node]==0) && (options[AZ_output] != AZ_none)) {
          AZ_printf_out("Warning:No previous definition for A_norm found. Was ");
          AZ_printf_out("AZ_iterate used\n\tpreviously and was the scaling object ");
          AZ_printf_out("passed in the same as for\n\tthis invokation of ");
          AZ_printf_out("AZ_iterate()?\n");
        }
      }


      if  ((options[AZ_ignore_scaling]) && (conv_info->scaling->A_norm != 0.0)){
        conv_info->A_norm = conv_info->scaling->A_norm;
      }
      else {
        if (Amat->matrix_type == AZ_MSR_MATRIX ||
            Amat->matrix_type == AZ_VBR_MATRIX) {
          conv_info->A_norm =
            AZ_gmax_matrix_norm(Amat->val, Amat->indx, Amat->bindx, Amat->rpntr,
                              Amat->cpntr, Amat->bpntr, proc_config,
                              Amat->data_org);
        }
        else {
          Amat->matrix_norm = Amat->matnorminf(Amat);
          if (Amat->matrix_norm <= 0.0) {
            if (proc_config[AZ_node]==0) {
              AZ_printf_out("Error: matrix is not MSR or VBR, and Amat->"
                "matrix_norm <= 0.0.\n");
            }
          }
          conv_info->A_norm = Amat->matrix_norm;
        }
      }

      if (fabs(conv_info->A_norm) < DBL_MIN) {
        AZ_p_error("Error: ||A||_infinity = 0\n\tSomething wrong with A\n",
                   proc_config[AZ_node]);
        exit(-1);
      }
    }

    else if (*r_avail) {
      dots[count] = DDOT_F77(&N, tr, &one, tr, &one);

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
          AZ_printf_out("Warning:No previous definition for A_norm found. Was ");
          AZ_printf_out("AZ_iterate used\n\tpreviously and was the scaling object ");
          AZ_printf_out("passed in the same as for\n\tthis invokation of ");
          AZ_printf_out("AZ_iterate()?\n");
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
                            AZ_ALLOC, AZ_SYS+az_iterate_id,
                            "temp in AZ_compute_global_scalars", &j);
    if (conv_info->not_initialized) {
      total_N = AZ_gsum_int(N, proc_config);
      conv_info->total_N = (double) total_N;
    }
    total_N = conv_info->total_N;

    for (i = 0; i < N; i++) temp[i] = tr[i] / w[i];

    dots[count] = DDOT_F77(&N, temp, &one, temp, &one);

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

  case AZTECOO_conv_test:

    if (v1 != NULL) {
      AZ_gdot_vec(1, dots, tmp, proc_config);
      *value = dots[0];
    }
    print_info = -1; /* Never print unless on PE 0 */
    if (proc_config[AZ_node] == 0) { 
      if (options[AZ_output]==AZ_none || options[AZ_output]==AZ_warnings)
	print_info = -1; /* Never print */
      else if (options[AZ_output]==AZ_last || options[AZ_output]==AZ_summary)
	print_info = 0; /* Print only if converged */
      else if (options[AZ_output]==AZ_all)
	print_info = 1; /* Print always */
      else if (conv_info->iteration%options[AZ_output]==0)
	print_info = 1; /* print this time */
      else 
	print_info = 0; /* Print only if converged */
    }
    conv_info->conv_function(conv_info->conv_object, conv_info->res_vec_object, 
			     conv_info->iteration,
                             r, print_info, conv_info->sol_updated,
                             &(conv_info->converged), &(conv_info->isnan), r_norm,
			     r_avail);

    *scaled_r_norm = *r_norm;

    break;

  default:
    if (proc_config[AZ_node] == 0) {
      (void) AZ_printf_err( "Error: Improper value, options[AZ_conv] = %d\n",
                     options[AZ_conv]);
    }
    exit(-1);
  }
  conv_info->not_initialized = AZ_FALSE;
  if (options[AZ_conv]!=AZTECOO_conv_test) {
    conv_info->converged = *scaled_r_norm < conv_info->epsilon;
    conv_info->isnan = (!(conv_info->converged) && !(*scaled_r_norm >= conv_info->epsilon));
  }

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
  int tmp_sol_updated;

  int r_avail = AZ_TRUE;

  /**************************** execution begins ******************************/

  /* calculate true residual */

  AZ_compute_residual(b, x, v, proc_config, Amat);

  /* compute scaled residual */

  tmp_sol_updated = conv_info->sol_updated;
  conv_info->sol_updated = 1; /* Solution is available here */
  AZ_compute_global_scalars(Amat, x, b, v, w,
                            actual_residual, scaled_r_norm, options, data_org,
                            proc_config,&r_avail, (double *) 0, (double *) 0,
                            (double *) 0, conv_info);

  conv_info->sol_updated = tmp_sol_updated;

} /* AZ_scale_true_residual */

struct AZ_CONVERGE_STRUCT *AZ_converge_create()
{
  struct AZ_CONVERGE_STRUCT *temp;

  temp = (struct AZ_CONVERGE_STRUCT *) AZ_allocate(sizeof(struct
                                                          AZ_CONVERGE_STRUCT));
  if (temp == NULL) {
    AZ_printf_err("AZ_converge_create: Not enough space\n");
    exit(1);
  }
  temp->r0_norm = 0.0;
  temp->A_norm  = 0.0;
  temp->b_norm  = 0.0;
  temp->scaling = NULL;
  temp->total_N       = 0;
  temp->not_initialized = AZ_TRUE;
  temp->isnan = 0;
  temp->converged = 0;

  return temp;
}
void AZ_converge_destroy(struct AZ_CONVERGE_STRUCT **temp)
{
  if (*temp != NULL) AZ_free(*temp);
  *temp = NULL;
}
