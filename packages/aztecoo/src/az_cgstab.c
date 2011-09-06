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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "az_aztec.h"
#include "az_blas_wrappers.h"

extern int az_iterate_id;

void AZ_pbicgstab(double b[], double x[], double weight[], int options[], 
	double params[],int proc_config[], double status[], AZ_MATRIX *Amat, 
	AZ_PRECOND *precond, struct AZ_CONVERGE_STRUCT *convergence_info)

/*******************************************************************************

  Vand der Vorst's (1990) variation of the Bi-Conjugate Gradient algorthm
  (Sonneveld (1984,1989)) to solve the nonsymmetric matrix problem Ax = b.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  b:               Right hand side of linear system.

  x:               On input, contains the initial guess. On output contains the
                   solution to the linear system.

  weight:          Vector of weights for convergence norm #4.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  status:          On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

  Amat:            Structure used to represent the matrix (see file az_aztec.h
                   and Aztec User's Guide).
  precond:         Structure used to represent the preconditionner
                   (see file az_aztec.h and Aztec User's Guide).

*******************************************************************************/

{

  /* local variables */

  register int    i;
  int             N, NN, one = 1, iter=1, r_avail = AZ_TRUE, j;
  int             precond_flag, print_freq, proc;
  int             brkdown_will_occur = AZ_FALSE;
  double          alpha = 1.0, beta, true_scaled_r=0.0;
  double          *v, *r, *rtilda, *p, *phat, *s, *shat;
  double          omega = 1.0, dot_vec[2], tmp[2], init_time = 0.0;
  double          rhonm1 = 1.0, rhon, sigma, brkdown_tol = DBL_EPSILON;
  double          scaled_r_norm= -1.0, actual_residual = -1.0, rec_residual= -1.0;
  double          dtemp;
  int          *data_org, str_leng, first_time = AZ_TRUE;
  char         label[64],suffix[32], prefix[64];


  /**************************** execution begins ******************************/

  sprintf(suffix," in cgstab%d",options[AZ_recursion_level]);
                                                           /* set string that will be used */
                                                           /* for manage_memory label      */
  /* set prefix for printing */

  str_leng = 0;
  for (i = 0; i < 16; i++) prefix[str_leng++] = ' ';
  for (i = 0 ; i < options[AZ_recursion_level]; i++ ) {
     prefix[str_leng++] = ' '; prefix[str_leng++] = ' '; prefix[str_leng++] = ' ';
     prefix[str_leng++] = ' '; prefix[str_leng++] = ' ';
  }
  prefix[str_leng] = '\0';

  data_org = Amat->data_org;

  /* pull needed values out of parameter arrays */

  N            = data_org[AZ_N_internal] + data_org[AZ_N_border];
  precond_flag = options[AZ_precond];
  proc         = proc_config[AZ_node];
  print_freq   = options[AZ_print_freq];

  /* Initialize some values in convergence info struct */
  convergence_info->print_info = print_freq;
  convergence_info->iteration = 0;
  convergence_info->sol_updated = 1; /* BiCGStab always updates solution */
  convergence_info->epsilon = params[AZ_tol]; /* Test against this */

  /* allocate memory for required vectors */

  NN     = N + data_org[AZ_N_external];
  if (NN == 0) NN++;   /* make sure everybody allocates something*/
  NN = NN + (NN%2);    /* make sure things are aligned for the   */
                       /* assembly coded matvec() on the Intel.  */

  sprintf(label,"phat%s",suffix);
  phat   = (double *) AZ_manage_memory(7*NN*sizeof(double), AZ_ALLOC, 
                                       AZ_SYS+az_iterate_id, label,&j);
  p      = &(phat[1*NN]);
  shat   = &(phat[2*NN]);   /* NOTE: phat and shat must be aligned */
                            /*       so that the assembly dgemv    */
                            /*       works on the paragon.         */
  s      = &(phat[3*NN]);
  r      = &(phat[4*NN]);
  rtilda = &(phat[5*NN]);
  v      = &(phat[6*NN]);

  AZ_compute_residual(b, x, r, proc_config, Amat);

  /* v, p <- 0 */

  for (i = 0; i < N; i++) v[i] = p[i] = 0.0;

  /* set rtilda */

  if (options[AZ_aux_vec] == AZ_resid)
    DCOPY_F77(&N, r, &one, rtilda, &one);
  else
    AZ_random_vector(rtilda, data_org, proc_config);

  /*
   * Compute a few global scalars:
   *     1) ||r||                corresponding to options[AZ_conv]
   *     2) scaled ||r||         corresponding to options[AZ_conv]
   *     3) rho = <rtilda, r>
   */

  AZ_compute_global_scalars(Amat, x, b, r,
                            weight, &rec_residual, &scaled_r_norm, options,
                            data_org, proc_config,&r_avail,r,rtilda, &rhon,
                            convergence_info);
  true_scaled_r = scaled_r_norm;

  if ((options[AZ_output] != AZ_none) && 
      (options[AZ_output] != AZ_last) &&
      (options[AZ_output] != AZ_warnings) &&
      (options[AZ_output] != AZ_summary) &&
      (options[AZ_conv]!=AZTECOO_conv_test) && (proc == 0))
    (void) AZ_printf_out("%siter:    0           residual = %e\n",prefix,scaled_r_norm);


  for (iter = 1; iter <= options[AZ_max_iter] && !(convergence_info->converged)
	 && !(convergence_info->isnan); iter++) {
    if (brkdown_will_occur) {
      AZ_scale_true_residual( x, b, v,
                             weight, &actual_residual, &true_scaled_r, options,
                             data_org, proc_config, Amat, convergence_info);

      AZ_terminate_status_print(AZ_breakdown, iter, status, rec_residual,
                                params, true_scaled_r, actual_residual, options,
                                proc_config);
      return;
    }

    beta = (rhon/rhonm1) * (alpha/omega);

    if (fabs(rhon) < brkdown_tol) { /* possible problem */
      if (AZ_breakdown_f(N, r, rtilda, rhon, proc_config))
        brkdown_will_occur = AZ_TRUE;
      else
        brkdown_tol = 0.1 * fabs(rhon);
    }

    rhonm1 = rhon;

    /* p    = r + beta*(p - omega*v)       */
    /* phat = M^-1 p                       */
    /* v    = A phat                       */

    dtemp = beta * omega;
    for (i = 0; i < N; i++) p[i] = r[i] + beta * p[i] - dtemp * v[i];
    DCOPY_F77(&N, p, &one, phat, &one);

    if (iter==1) init_time = AZ_second();
    if (precond_flag)
      precond->prec_function(phat,options,proc_config,params,Amat,precond);

    if (iter==1) status[AZ_first_precond] = AZ_second() - init_time;

    Amat->matvec(phat, v, Amat, proc_config);

    sigma = AZ_gdot(N, rtilda, v, proc_config);

    if (fabs(sigma) < brkdown_tol) { /* possible problem */
      if (AZ_breakdown_f(N, rtilda, v, sigma, proc_config)) {

        /* break down */

        AZ_scale_true_residual( x, b, v,
                               weight, &actual_residual, &true_scaled_r,
                               options, data_org,proc_config, Amat,
			       convergence_info);

        AZ_terminate_status_print(AZ_breakdown, iter, status, rec_residual,
                                  params, true_scaled_r, actual_residual,
                                  options, proc_config);
        return;
      }
      else brkdown_tol = 0.1 * fabs(sigma);
    }

    alpha = rhon / sigma;

    /* s = r - alpha*v                     */
    /* shat = M^-1 s                       */
    /* r = A shat (r is a tmp here for t ) */

    for (i = 0; i < N; i++) s[i] = r[i] - alpha * v[i];
    DCOPY_F77(&N, s, &one, shat, &one);

    if (precond_flag)    
      precond->prec_function(shat,options,proc_config,params,Amat,precond);

    Amat->matvec(shat, r, Amat, proc_config);


    /* omega = (t,s)/(t,t) with r = t */

    dot_vec[0] = DDOT_F77(&N, r, &one, s, &one);
    dot_vec[1] = DDOT_F77(&N, r, &one, r, &one);
    AZ_gdot_vec(2, dot_vec, tmp, proc_config);

    if (fabs(dot_vec[1]) < DBL_MIN) {
      omega = 0.0;
      brkdown_will_occur = AZ_TRUE;
    }
    else omega = dot_vec[0] / dot_vec[1];

    /* x = x + alpha*phat + omega*shat */
    /* r = s - omega*r */

    DAXPY_F77(&N, &alpha, phat, &one, x, &one);
    DAXPY_F77(&N, &omega, shat, &one, x, &one);

    for (i = 0; i < N; i++) r[i] = s[i] - omega * r[i];

    /*
     * Compute a few global scalars:
     *     1) ||r||                corresponding to options[AZ_conv]
     *     2) scaled ||r||         corresponding to options[AZ_conv]
     *     3) rho = <rtilda, r>
     */

    AZ_compute_global_scalars(Amat, x, b, r,
                              weight, &rec_residual, &scaled_r_norm, options,
                              data_org, proc_config, &r_avail, r, rtilda, &rhon,
                              convergence_info);

    if ( (iter%print_freq == 0) && proc == 0)
      (void) AZ_printf_out("%siter: %4d           residual = %e\n",prefix,iter,
                     scaled_r_norm);

    /* convergence tests */

    if (options[AZ_check_update_size] & convergence_info->converged) {
      dtemp = alpha/omega;
      DAXPY_F77(&N, &dtemp, phat, &one, shat, &one);
      convergence_info->converged = AZ_compare_update_vs_soln(N, -1.,omega, shat, x,
                                         params[AZ_update_reduction],
                                         options[AZ_output], proc_config, &first_time);
    }

    if (convergence_info->converged) {
      AZ_scale_true_residual(x, b, v,
                             weight, &actual_residual, &true_scaled_r, options,
                             data_org, proc_config, Amat, convergence_info);


      /*
       * Note: epsilon and params[AZ_tol] may not be equal due to a previous
       *       call to AZ_get_new_eps().
       */

      if (!(convergence_info->converged) && options[AZ_conv]!=AZTECOO_conv_test) {
	if (AZ_get_new_eps(&convergence_info->epsilon, scaled_r_norm, true_scaled_r,
			   options, proc_config) == AZ_QUIT) {

	  /*
	   * Computed residual has converged, actual residual has not converged,
	   * AZ_get_new_eps() has decided that it is time to quit.
	   */
	  
	  AZ_terminate_status_print(AZ_loss, iter, status, rec_residual, params,
				    true_scaled_r, actual_residual, options,
				    proc_config);
	  return;
	}
      }
    }
  }

  iter--;
  if ( (iter%print_freq != 0) && (proc == 0) && (options[AZ_output] != AZ_none)
       && (options[AZ_output] != AZ_warnings) &&
       (options[AZ_conv]!=AZTECOO_conv_test))
    (void) AZ_printf_out("%siter: %4d           residual = %e\n", prefix,iter,
                   scaled_r_norm);

  /* check if we exceeded maximum number of iterations */

  if (convergence_info->converged) {
    i = AZ_normal;
    scaled_r_norm = true_scaled_r;
  }
  else if (convergence_info->isnan) i = AZ_breakdown;
  else
    i = AZ_maxits;

  AZ_terminate_status_print(i, iter, status, rec_residual, params,
                            scaled_r_norm, actual_residual, options,
                            proc_config);

} /* bicgstab */
