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
 * $RCSfile: az_gmresr.c,v $
 *
 * $Author: william $
 *
 * $Date: 2006/07/20 23:08:48 $
 *
 * $Revision: 1.11 $
 *
 * $Name:  $
 *====================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "az_aztec.h"
#include "az_blas_wrappers.h"

extern int az_iterate_id;

void AZ_pgmresr(double b[], double x[],double weight[], int options[],
	double params[], int proc_config[], double status[], AZ_MATRIX *Amat, 
	AZ_PRECOND *precond, struct AZ_CONVERGE_STRUCT *convergence_info)

/*******************************************************************************

  This routine uses Saad's restarted Genralized Minimum Residual method to solve
  the nonsymmetric matrix problem Ax = b.

  IMPORTANT NOTE: While the 2-norm of the gmres residual is available, the
  actual residual is not normally computed as part of the gmres algorithm. Thus,
  if the user uses a convergence condition (see AZ_gmres_global_scalars()) that
  is based on the 2-norm of the residual there is no need to compute the
  residual (i.e. r_avail = AZ_FALSE). However, if another norm of r is
  requested, AZ_gmres_global_scalars() sets r_avail = AZ_TRUE and the algorithm
  computes the residual.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  Amat:            Structure used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  b:               Right hand side of linear system.

  x:               On input, contains the initial guess. On output contains the
                   solution to the linear system.

  weight:          Vector of weights for convergence norm #4.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file Aztec User's Guide).

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
*******************************************************************************/

{

  /* local variables */

  register int k;
  int          i, N, NN, converged, one = 1, iter, r_avail = AZ_FALSE;
  int          print_freq, proc, kspace;
  double     **UU, **CC, *dots, *tmp, *res;
  double       dble_tmp, r_2norm = 1.0, epsilon;
  double       rec_residual, scaled_r_norm, true_scaled_r=0.0;
  double       actual_residual = -1.0, minus_alpha, alpha;
  double       *dummy = (double *) 0;
  double       *UUblock, *CCblock;
  int          mm, ii;
  char         label[64],suffix[32], prefix[64];
  int          *data_org, str_leng, first_time = AZ_TRUE;
  double       doubleone = 1.0, minusone = -1.0, init_time = 0.0;
char *T = "T";
char *T2 = "N";


  /**************************** execution begins ******************************/

  sprintf(suffix," in gmresr%d",options[AZ_recursion_level]);
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
  epsilon      = params[AZ_tol];
  proc         = proc_config[AZ_node];
  print_freq   = options[AZ_print_freq];
  kspace       = options[AZ_kspace];

  /* Initialize some values in convergence info struct */
  convergence_info->print_info = print_freq;
  convergence_info->iteration = 0;
  convergence_info->sol_updated = 0; /* GMRES seldom updates solution */
  convergence_info->epsilon = params[AZ_tol];

  /* allocate memory for required vectors */

  NN    = kspace  + 1;
  /* +1: make sure everybody allocates something */

  sprintf(label,"dots%s",suffix);
  dots  = AZ_manage_memory(2*NN*sizeof(double), AZ_ALLOC,AZ_SYS+az_iterate_id,label,&i);
  tmp   = &(dots[NN]);
  sprintf(label,"CC%s",suffix);
  CC    = (double **) AZ_manage_memory(2*NN*sizeof(double *),
                                       AZ_ALLOC,AZ_SYS+az_iterate_id,label,&i);
  UU    = &(CC[NN]);

  NN    = N + data_org[AZ_N_external];
  if (NN == 0) NN++; /* make sure everybody allocates something */
  NN = NN + (NN%2);  /* make sure things are aligned for intel  */


  sprintf(label,"UUblock%s",suffix);
  UUblock = AZ_manage_memory(2*NN*kspace*sizeof(double),
                             AZ_ALLOC, AZ_SYS+az_iterate_id,label, &i);
  for (k = 0; k < kspace; k++) UU[k] = &(UUblock[k*NN]);
  CCblock = &(UUblock[kspace*NN]);
  for (k = 0; k < kspace; k++) CC[k] = &(CCblock[k*NN]);

  sprintf(label,"res%s",suffix);
  res = AZ_manage_memory(NN*sizeof(double),AZ_ALLOC,AZ_SYS+az_iterate_id,label,&i);

  AZ_compute_residual(b, x, res, proc_config, Amat);

  /*
   * Compute a few global scalars:
   *     1) ||r||                corresponding to options[AZ_conv]
   *     2) scaled ||r||         corresponding to options[AZ_conv]
   */
  r_2norm = DDOT_F77(&N, res, &one, res, &one);
  AZ_gdot_vec(1, &r_2norm, &rec_residual, proc_config);  
  r_2norm = sqrt(r_2norm);
  rec_residual = r_2norm;

  AZ_compute_global_scalars(Amat, x, b, res,
                          weight, &rec_residual, &scaled_r_norm, options,
                          data_org, proc_config, &r_avail, NULL, NULL, NULL,
                          convergence_info);
  r_2norm = rec_residual;

  converged = scaled_r_norm < epsilon;

  if ( (options[AZ_output] != AZ_none) && 
       (options[AZ_output] != AZ_last) &&
       (options[AZ_output] != AZ_summary) &&
       (options[AZ_output] != AZ_warnings) && (proc == 0) )
    (void) AZ_printf_out("%siter:    0           residual = %e\n",
                           prefix,scaled_r_norm);

  iter = 0;
/*rst change  while (!converged && iter < options[AZ_max_iter]) { */
  while (!(convergence_info->converged) && iter < options[AZ_max_iter] && !(convergence_info->isnan)) {

    convergence_info->iteration = iter;
    i = 0;

/*rst change   while (i < kspace && !converged && iter < options[AZ_max_iter]) { */
    while (i < kspace && !(convergence_info->converged) && iter < options[AZ_max_iter]
           && !(convergence_info->isnan)) {

      iter++;
    convergence_info->iteration = iter;


      /* v_i+1 = A M^-1 v_i */

      DCOPY_F77(&N, res , &one, UU[i], &one);

      if (iter == 1) init_time = AZ_second();

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
        /* Start timer. */
      static int precID = -1;
      precID = Teuchos_startTimer( "AztecOO: Operation Prec*x", precID );
#endif
#endif
      precond->prec_function(UU[i],options,proc_config,params,Amat,precond);
#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Stop timer. */
      Teuchos_stopTimer( precID );
#endif
#endif
      if (iter == 1) status[AZ_first_precond] = AZ_second() - init_time;

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Start timer. */
      static int matvecID = -1;
      matvecID = Teuchos_startTimer( "AztecOO: Operation Op*x", matvecID );
#endif
#endif
      Amat->matvec(UU[i], CC[i], Amat, proc_config);
#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Stop timer. */
      Teuchos_stopTimer( matvecID );
#endif
#endif

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Start the timer. */
      static int orthoID = -1;
      orthoID = Teuchos_startTimer( "AztecOO: Orthogonalization", orthoID );
#endif
#endif

      /* Gram-Schmidt orthogonalization */

      if (!options[AZ_orthog]) { /* classical  (stabilized) */
         for (ii = 0 ; ii < 2 ; ii++ ) {
            dble_tmp = 0.0; mm = i;
            if (N == 0) for (k = 0 ; k < i ; k++) dots[k] = 0.0;
#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Start the timer. */
      static int orthoInnerProdID = -1;
      orthoInnerProdID = Teuchos_startTimer( "AztecOO: Ortho (Inner Product)", orthoInnerProdID );
#endif
#endif
            DGEMV_F77(CHAR_MACRO(T[0]), &N, &mm, &doubleone, CCblock, &NN, CC[i], 
                   &one, &dble_tmp, dots, &one);

            AZ_gdot_vec(i, dots, tmp, proc_config);
#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      Teuchos_stopTimer( orthoInnerProdID );
#endif
#endif

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Start the timer. */
      static int orthoUpdateID = -1;
      orthoUpdateID = Teuchos_startTimer( "AztecOO: Ortho (Update)", orthoUpdateID );
#endif
#endif
            DGEMV_F77(CHAR_MACRO(T2[0]), &N, &mm, &minusone, CCblock, &NN, dots, 
                   &one, &doubleone, CC[i], &one);
            DGEMV_F77(CHAR_MACRO(T2[0]), &N, &mm, &minusone, UUblock, &NN, dots,
                   &one, &doubleone, UU[i], &one);
#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      Teuchos_stopTimer( orthoUpdateID );
#endif
#endif
         }
      }
      else {                    /* modified */
        for (k = 0; k < i; k++) {
          alpha = AZ_gdot(N, CC[k], CC[i], proc_config);
          minus_alpha = -alpha;
          DAXPY_F77(&N, &minus_alpha, CC[k], &one, CC[i], &one);
          DAXPY_F77(&N, &minus_alpha, UU[k], &one, UU[i], &one);
        }
      }

      /* normalize vector */

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      static int orthoNormID = -1;
      orthoNormID = Teuchos_startTimer( "AztecOO: Ortho (Norm)", orthoNormID );
#endif
#endif
      dble_tmp = sqrt(AZ_gdot(N, CC[i], CC[i], proc_config));
#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      Teuchos_stopTimer( orthoNormID );
#endif
#endif

      if (dble_tmp  > DBL_EPSILON*r_2norm)
        dble_tmp  = 1.0 / dble_tmp;
      else
        dble_tmp = 0.0;

      DSCAL_F77(&N, &dble_tmp, CC[i], &one);
      DSCAL_F77(&N, &dble_tmp, UU[i], &one);

      dble_tmp = AZ_gdot(N, CC[i], res, proc_config);
      DAXPY_F77(&N, &dble_tmp, UU[i], &one, x, &one);
      dble_tmp = -dble_tmp;
      DAXPY_F77(&N, &dble_tmp, CC[i], &one, res, &one);

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Stop the timer. */
      Teuchos_stopTimer( orthoID );
#endif
#endif

      /* determine residual norm & test convergence */

      r_2norm      = sqrt(AZ_gdot(N, res, res, proc_config));
      rec_residual = r_2norm;

      /*
       * Compute a few global scalars:
       *     1) ||r||                corresponding to options[AZ_conv]
       *     2) scaled ||r||         corresponding to options[AZ_conv]
       * NOTE: if r_avail = AZ_TRUE or AZ_FIRST is passed in, we perform
       * step 1), otherwise ||r|| is taken as rec_residual.
       */

      AZ_compute_global_scalars(Amat, x, b, res,
                              weight, &rec_residual, &scaled_r_norm, options,
                              data_org, proc_config, &r_avail, dummy, dummy,
                              dummy, convergence_info);

      converged = scaled_r_norm < epsilon;

/*rst change      if ( (iter%print_freq == 0) && proc == 0) */
      if ( (iter%print_freq == 0) &&
           (options[AZ_conv]!=AZTECOO_conv_test) && proc == 0)
        (void) AZ_printf_out("%siter: %4d           residual = %e\n",prefix,iter,
                       scaled_r_norm);

      i++;      /* subspace dim. counter dim(K) = i - 1 */
#ifdef out
      if (options[AZ_check_update_size] & converged)
         converged = AZ_compare_update_vs_soln(N, -1.,dble_tmp, UU[i-1], x,
                                           params[AZ_update_reduction],
                                           options[AZ_output], proc_config, &first_time);



      if (converged) {

        /* compute true residual using 'v[kspace]' as a temporary vector */

        AZ_scale_true_residual(x, b,
                               res, weight, &actual_residual, &true_scaled_r,
                               options, data_org, proc_config, Amat,
			       convergence_info);

        converged = true_scaled_r < params[AZ_tol];

        if (!converged && (AZ_get_new_eps(&epsilon, scaled_r_norm,
                                          true_scaled_r,
                                          options, proc_config) == AZ_QUIT)) {

          /*
           * Computed residual has converged, actual residual has not
           * converged, AZ_get_new_eps() has decided that it is time to quit.
           */

          AZ_terminate_status_print(AZ_loss, iter, status, rec_residual, params,
                                    true_scaled_r, actual_residual, options,
                                    proc_config);
          return;
        }
      }
#endif
    }
  }

  if ( (iter%print_freq != 0) && (proc == 0) && (options[AZ_output] != AZ_none)
       && (options[AZ_output] != AZ_warnings))
    (void) AZ_printf_out("%siter: %4d           residual = %e\n",
		   prefix,iter, scaled_r_norm);


  if (convergence_info->converged) {
    i = AZ_normal;
    scaled_r_norm = true_scaled_r;
  }
  else if (convergence_info->isnan) i = AZ_breakdown;
  else i = AZ_maxits;

  AZ_terminate_status_print(i, iter, status, rec_residual, params,
                            scaled_r_norm, actual_residual, options,
                            proc_config);

#ifdef out
  /* check if we exceeded maximum number of iterations */

  if (converged) {
    i = AZ_normal;
    scaled_r_norm = true_scaled_r;
  }
  else
    i = AZ_maxits;

  AZ_terminate_status_print(i, iter, status, rec_residual, params,
                            scaled_r_norm, actual_residual, options,
                            proc_config);

#endif
} /* AZ_pgmres */

