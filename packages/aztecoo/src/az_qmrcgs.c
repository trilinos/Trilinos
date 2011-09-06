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

void AZ_pqmrs(double b[], double x[], double weight[], int options[], 
	double params[], int proc_config[], double status[], AZ_MATRIX *Amat, 
	AZ_PRECOND *precond, struct AZ_CONVERGE_STRUCT *convergence_info)

/*******************************************************************************

  Freund's transpose free QMR routine to solve the nonsymmetric matrix problem
  Ax = b. NOTE: this routine differs from Freund's paper in that we compute
  ubar (= M^-1 u ) and qbar (= M^-1 q) instead of u and q defined in Freund's
  paper.

  IMPORTANT NOTE: While an estimate of the 2-norm of the qmr residual is
  available (see comment below), the actual qmr residual is not normally
  computed as part of the qmr algorithm. Thus, if the user uses a convergence
  condition (see AZ_compute_global_scalars()) that is based on the 2-norm of the
  residual there is no need to compute the residual (i.e. r_avail = AZ_FALSE).
  However, if another norm of r is requested, AZ_compute_global_scalars() will
  set r_avail = AZ_TRUE and the algorithm will compute the residual.

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

  Amat:            Structure used to represent the matrix (see az_aztec.h
                   and Aztec User's Guide).
  Oprecond:         Structure used to represent the preconditionner
                   (see file az_aztec.h and Aztec User's Guide).

*******************************************************************************/

{

  /* local variables */

  register int    i;
  int             N, NN, one = 1, iter= 1,r_avail = AZ_FALSE, j;
  int             precond_flag, print_freq, proc;
  int             brkdown_will_occur = AZ_FALSE;
  double          alpha, beta = 0.0, true_scaled_r=0.0;
  double          *ubar, *v, *r_cgs, *rtilda, *Aubar, *qbar, *Aqbar, *d, *Ad = NULL;
  double          rhonm1, rhon, est_residual, actual_residual = -1.0;
  double          scaled_r_norm, sigma, brkdown_tol = DBL_EPSILON;
  double          omega, c, norm_r_n_cgs, norm_r_nm1_cgs;
  double          tau_m, nu_m, eta_m, init_time = 0.0;
  double          tau_mm1, nu_mm1 = 0.0, eta_mm1 = 0.0, doubleone = 1.0;
  register double dtemp;
  double          W_norm = 0.0;
  int             offset = 0;
  int          *data_org, str_leng, first_time = AZ_TRUE;
  char         label[64],suffix[32], prefix[64];

  /**************************** execution begins ******************************/

  sprintf(suffix," in qmrcgs%d",options[AZ_recursion_level]);
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
  convergence_info->sol_updated = 1; /* QMRCGS always updates solution */
  convergence_info->epsilon = params[AZ_tol]; /* Test against this */

  /* allocate memory for required vectors */

  NN     = N + data_org[AZ_N_external];
  if (NN == 0) NN++; /* make sure everyone allocates something */
  NN = NN + (NN%2);
      /* make sure things are aligned on double words for paragon */


  sprintf(label,"ubar%s",suffix);
  ubar   = (double *) AZ_manage_memory(8*NN*sizeof(double),
				      AZ_ALLOC,AZ_SYS+az_iterate_id,label,&j);
  v      = &(ubar[1*NN]);
  Aubar  = &(ubar[2*NN]);
  d      = &(ubar[3*NN]);
  qbar   = &(ubar[4*NN]);
  rtilda = &(ubar[5*NN]);
  Aqbar  = &(ubar[6*NN]);
  r_cgs  = &(ubar[7*NN]);

  AZ_compute_residual(b, x, r_cgs, proc_config, Amat);

  /* d, qbar, Aqbar, v = 0 */

  for (i = 0; i < N; i++)
    d[i] = qbar[i] = Aqbar[i] = v[i] = 0.0;

  /* set rtilda */

  if (options[AZ_aux_vec] == AZ_resid) DCOPY_F77(&N, r_cgs, &one, rtilda, &one);
  else AZ_random_vector(rtilda, data_org, proc_config);

  /*
   * Compute a few global scalars:
   *     1) ||r_cgs||              corresponding to options[AZ_conv]
   *     2) scaled ||r_cgs||       corresponding to options[AZ_conv]
   *     3) rhon = <rtilda, r_cgs>
   * Note: step 1) is performed if r_avail = AZ_TRUE on entry or
   *       AZ_FIRST_TIME is passed in. Otherwise, ||r_cgs|| is taken as
   *       est_residual.
   */

  /* Change to support AztecOO_StatusTest:
     Even though AZ_compute_global_scalars can compute est_residual for us,
     we want to compute the initial residual prior to calling 
     this function because we need this value passed in to 
     AZ_compute_global_scalars in order to be consistent with subsequent calls
     to the same function later.  
  */
  scaled_r_norm = DDOT_F77(&N, r_cgs, &one, r_cgs, &one);
  AZ_gdot_vec(1, &scaled_r_norm, &est_residual, proc_config); 
  scaled_r_norm = sqrt(scaled_r_norm);
  est_residual = scaled_r_norm;

  AZ_compute_global_scalars(Amat, x, b, r_cgs,
                            weight, &est_residual, &scaled_r_norm, options,
                            data_org, proc_config, &r_avail, r_cgs, rtilda,
                            &rhon, convergence_info);
  true_scaled_r = scaled_r_norm;

  if ((options[AZ_output] != AZ_none) && 
      (options[AZ_output] != AZ_last) &&
      (options[AZ_output] != AZ_warnings) &&
      (options[AZ_output] != AZ_summary) &&
      (options[AZ_conv]!=AZTECOO_conv_test) && (proc == 0))
    (void) AZ_printf_out("%siter:    0           residual = %e\n",prefix,scaled_r_norm);

  norm_r_nm1_cgs = est_residual;
  tau_mm1        = norm_r_nm1_cgs;
  rhonm1         = rhon;

  /* Set up aux-vector if we need to compute the qmr residual */

  /* We always want the residual for AztecOO tests */
  if (r_avail || (options[AZ_conv]==AZTECOO_conv_test)) {
    sprintf(label,"Ad%s",suffix);
    Ad = (double *) AZ_manage_memory(NN*sizeof(double),AZ_ALLOC,
				     AZ_SYS+az_iterate_id, label, &j);
    for (i = 0; i < N; i++) Ad[i] = 0.0;
  }


  for (iter = 1; iter <= options[AZ_max_iter] && !(convergence_info->converged) && 
	 !(convergence_info->isnan); iter++) {
    convergence_info->iteration = iter;
    
    if (fabs(rhon) < brkdown_tol) { /* possible breakdown problem */

      if (AZ_breakdown_f(N, r_cgs, rtilda, rhon, proc_config))
        brkdown_will_occur = AZ_TRUE;
      else
        brkdown_tol = 0.1 * fabs(rhon);
    }

    /* ubar  = M^-1 r_cgs + beta*qbar               */
    /* Aubar = A ubar                               */
    /* v     = A ubar + beta ( A qbar + beta pnm1 ) */
    /*       = Aubar  + beta ( Aqbar +  beta v)     */

    DCOPY_F77(&N, r_cgs, &one, ubar, &one);

    if (iter==1) init_time = AZ_second();

    if (precond_flag)
    precond->prec_function(ubar,options,proc_config,params,Amat,precond);

    if (iter==1) status[AZ_first_precond] = AZ_second() - init_time;

    for (i = 0; i < N; i++) ubar[i] = ubar[i] + beta * qbar[i];

    Amat->matvec(ubar, Aubar, Amat, proc_config);

    DAXPY_F77(&N, &beta, v, &one, Aqbar, &one);
    for (i = 0; i < N; i++) v[i] = Aubar[i] + beta * Aqbar[i];

    sigma = AZ_gdot(N, rtilda, v, proc_config);

    if (fabs(sigma) < brkdown_tol) { /* possible problem */
      if (AZ_breakdown_f(N, rtilda, v, sigma, proc_config)) {

        /* break down */

        AZ_scale_true_residual(x, b, v,
                               weight, &actual_residual, &true_scaled_r,
                               options, data_org, proc_config, Amat,
			       convergence_info);

        AZ_terminate_status_print(AZ_breakdown, iter, status, est_residual,
                                  params, true_scaled_r, actual_residual,
                                  options, proc_config);
        return;
      }
      else brkdown_tol = 0.1 * fabs(sigma);
    }

    alpha = rhon / sigma;

    /* qbar  = ubar - alpha* M^-1 v            */
    /* Aqbar = A qbar                          */
    /* r_cgs = r_cgs - alpha (A ubar + A qbar) */
    /*       = r_cgs - alpha (Aubar + Aqbar)   */

    DCOPY_F77(&N, v, &one, qbar, &one);

    if (precond_flag)
    precond->prec_function(qbar,options,proc_config,params,Amat,precond);

    for (i = 0; i < N; i++) qbar[i] = ubar[i] - alpha * qbar[i];
    Amat->matvec(qbar, Aqbar, Amat, proc_config);

    for (i = 0; i < N; i++) r_cgs[i] = r_cgs[i] - alpha*(Aubar[i] + Aqbar[i]);

    /* QMRS scaling and iterates weights 5.11 */

    norm_r_n_cgs = sqrt(AZ_gdot(N, r_cgs, r_cgs, proc_config));

    /* m is odd in Freund's paper */

    omega = sqrt(norm_r_nm1_cgs * norm_r_n_cgs);
    nu_m  = omega / tau_mm1;
    c     = 1.0 / sqrt(1.0 + nu_m * nu_m);
    tau_m = tau_mm1 * nu_m * c;
    eta_m = c * c * alpha;

    if (brkdown_will_occur) {
      AZ_scale_true_residual(x, b, v,
                             weight, &actual_residual, &true_scaled_r, options,
                             data_org, proc_config, Amat, convergence_info);

      AZ_terminate_status_print(AZ_breakdown, iter, status, est_residual,
                                params, true_scaled_r, actual_residual, options,
                                proc_config);
      return;
    }

    dtemp = nu_mm1 *nu_mm1 * eta_mm1 / alpha;
    for (i = 0; i < N; i++) d[i] = ubar[i] + dtemp * d[i];
    DAXPY_F77(&N, &eta_m, d, &one, x, &one); /* x = x - eta_m d  */

    if (r_avail || (options[AZ_conv]==AZTECOO_conv_test)) {
      for (i = 0; i < N; i++) Ad[i] = Aubar[i] + dtemp * Ad[i];
    }

    /* save some values */

    eta_mm1 = eta_m;  tau_mm1        = tau_m;
    nu_mm1  = nu_m;   norm_r_nm1_cgs = norm_r_n_cgs;

    /* m is even in Freund's paper */

    omega = norm_r_n_cgs;

    if (tau_mm1 == 0.0) nu_m = 0.0;
    else                nu_m  = omega / tau_mm1;

    c     = 1.0 / sqrt(1.0 + nu_m * nu_m);
    tau_m = tau_mm1 * nu_m * c;

    if (options[AZ_check_update_size]) {
       eta_m = eta_m/(c*c*alpha);
       for (i = 0; i < N; i++) ubar[i] = eta_m*d[i];
    }
    eta_m = c * c * alpha;

    dtemp = nu_mm1 * nu_mm1 * eta_mm1 / alpha;
    for (i = 0; i < N; i++) d[i] = qbar[i] + dtemp * d[i];
    DAXPY_F77(&N, &eta_m, d, &one, x, &one); /* x = x - eta_m d  */

    if (r_avail || (options[AZ_conv]==AZTECOO_conv_test)) {
      for (i = 0; i < N; i++) Ad[i] = Aqbar[i] + dtemp * Ad[i];
    }

    /* save some values */

    eta_mm1 = eta_m;  tau_mm1        = tau_m;
    nu_mm1  = nu_m;   norm_r_nm1_cgs = norm_r_n_cgs;
    rhonm1  = rhon;

    if (r_avail || (options[AZ_conv]==AZTECOO_conv_test)) {
      for (i = 0; i < N; i++) Aubar[i] = r_cgs[i] - (eta_m - alpha) * Ad[i];
      if (options[AZ_conv]==AZTECOO_conv_test) {
	
	scaled_r_norm = DDOT_F77(&N, Aubar, &one, Aubar, &one);
	AZ_gdot_vec(1, &scaled_r_norm, &est_residual, proc_config);  
	est_residual = scaled_r_norm;
      }
      
      /* Note: Aubar temporarily holds qmr residual */
      
    }
    else {

      /*
       * We want to estimate the 2-norm of the qmr residual. Freund gives the
       * bound ||r|| <= tau_m * sqrt(2*iter+1).  We use this bound until we get
       * close to the solution. At that point we compute the real residual norm
       * and use this to estimate the norm of ||W|| in Freund's paper.
       */

      dtemp = sqrt((double) (2 * iter + 1));

      if ((scaled_r_norm < convergence_info->epsilon * dtemp) && !offset) {
        AZ_scale_true_residual(x, b,
                               Aubar, weight, &actual_residual, &true_scaled_r,
                               options, data_org, proc_config, Amat,
			       convergence_info);

        if (tau_m != 0.0) W_norm = actual_residual / tau_m;
        if (W_norm < 1.0) W_norm = 1.0;

        offset       = 2 * iter + 1;
        est_residual = actual_residual;
      }
      else est_residual = sqrt((double)(2 * iter + 1 - offset) +
                               W_norm * W_norm) * tau_m;
    }

    /*
     * Compute a few global scalars:
     *     1) ||r||                corresponding to options[AZ_conv]
     *     2) scaled ||r||         corresponding to options[AZ_conv]
     *     3) rhon = <rtilda, r_cgs>
     * Note: step 1) is performed if r_avail = AZ_TRUE or AZ_FIRST_TIME
     *       is passed in. Otherwise, ||r|| is taken as est_residual.
     */

    AZ_compute_global_scalars(Amat, x, b,
                              Aubar, weight, &est_residual, &scaled_r_norm,
                              options, data_org, proc_config, &r_avail, rtilda,
                              r_cgs, &rhon, convergence_info);

    if ( (iter%print_freq == 0)  &&
	 (options[AZ_conv]!=AZTECOO_conv_test) && proc == 0 )
      (void) AZ_printf_out("%siter: %4d           residual = %e\n",prefix,iter,
                     scaled_r_norm);

    /* convergence tests */

    if (options[AZ_check_update_size] & convergence_info->converged) {
      DAXPY_F77(&N, &doubleone , d, &one, ubar, &one); 
      convergence_info->converged = AZ_compare_update_vs_soln(N, -1.,eta_m, ubar, x,
                                           params[AZ_update_reduction],
                                           options[AZ_output], proc_config, &first_time);
    }

    if (convergence_info->converged) {
      AZ_scale_true_residual(x, b, Aubar,
                             weight, &actual_residual, &true_scaled_r, options,
                             data_org, proc_config, Amat,convergence_info);


      /*
       * Note: epsilon and params[AZ_tol] may not be equal due to a previous
       * call to AZ_get_new_eps().
       */

      if (!(convergence_info->converged)  && options[AZ_conv]!=AZTECOO_conv_test) {
	
	if (AZ_get_new_eps(&convergence_info->epsilon, scaled_r_norm, true_scaled_r,
			   options, proc_config) == AZ_QUIT) {
	  
	  /*
	   * Computed residual has converged, actual residual has not converged,
	   * AZ_get_new_eps() has decided that it is time to quit.
	   */
	  
	  AZ_terminate_status_print(AZ_loss, iter, status, est_residual, params,
				    true_scaled_r, actual_residual, options,
				    proc_config);
	  return;
	}
      }
    }
    beta = rhon / rhonm1;
  }

  iter--;

  if ( (iter%print_freq != 0) &&
       (options[AZ_conv]!=AZTECOO_conv_test)  && (proc == 0) && (options[AZ_output] != AZ_none)
       && (options[AZ_output] != AZ_warnings))
    (void) AZ_printf_out("%siter: %4d           residual = %e\n",prefix,iter,
                   scaled_r_norm);

  /* check if we exceeded maximum number of iterations */

  if (convergence_info->converged) {
    i = AZ_normal;
    scaled_r_norm = true_scaled_r;
  }
  else if (convergence_info->isnan) i = AZ_breakdown;
  else
    i = AZ_maxits;

  AZ_terminate_status_print(i, iter, status, est_residual, params,
                            scaled_r_norm, actual_residual, options,
                            proc_config);

} /* pqmrs */
