/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
  dots  = AZ_manage_memory(2*NN*sizeof(double), AZ_ALLOC, AZ_SYS,label,&i);
  tmp   = &(dots[NN]);
  sprintf(label,"CC%s",suffix);
  CC    = (double **) AZ_manage_memory(2*NN*sizeof(double *),
                                       AZ_ALLOC,AZ_SYS,label,&i);
  UU    = &(CC[NN]);

  NN    = N + data_org[AZ_N_external];
  if (NN == 0) NN++; /* make sure everybody allocates something */
  NN = NN + (NN%2);  /* make sure things are aligned for intel  */


  sprintf(label,"UUblock%s",suffix);
  UUblock = AZ_manage_memory(2*NN*kspace*sizeof(double),
                             AZ_ALLOC, AZ_SYS,label, &i);
  for (k = 0; k < kspace; k++) UU[k] = &(UUblock[k*NN]);
  CCblock = &(UUblock[kspace*NN]);
  for (k = 0; k < kspace; k++) CC[k] = &(CCblock[k*NN]);

  sprintf(label,"res%s",suffix);
  res = AZ_manage_memory(NN*sizeof(double),AZ_ALLOC,AZ_SYS,label,&i);

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

  if ( (options[AZ_output] != AZ_none) && (options[AZ_output] != AZ_last) &&
       (options[AZ_output] != AZ_warnings) && (proc == 0) )
    (void) fprintf(stdout, "%siter:    0           residual = %e\n",
                           prefix,scaled_r_norm);

  iter = 0;
  while (!converged && iter < options[AZ_max_iter]) {

    i = 0;

    while (i < kspace && !converged && iter < options[AZ_max_iter]) {
      iter++;

      /* v_i+1 = A M^-1 v_i */

      DCOPY_F77(&N, res , &one, UU[i], &one);

      if (iter == 1) init_time = AZ_second();

      precond->prec_function(UU[i],options,proc_config,params,Amat,precond);
      if (iter == 1) status[AZ_first_precond] = AZ_second() - init_time;
      Amat->matvec(UU[i], CC[i], Amat, proc_config);


      /* Gram-Schmidt orthogonalization */

      if (!options[AZ_orthog]) { /* classical  (stabilized) */
         for (ii = 0 ; ii < 2 ; ii++ ) {
            dble_tmp = 0.0; mm = i;
            if (N == 0) for (k = 0 ; k < i ; k++) dots[k] = 0.0;
            DGEMV_F77(CHAR_MACRO(T[0]), &N, &mm, &doubleone, CCblock, &NN, CC[i], 
                   &one, &dble_tmp, dots, &one);

            AZ_gdot_vec(i, dots, tmp, proc_config);

            DGEMV_F77(CHAR_MACRO(T2[0]), &N, &mm, &minusone, CCblock, &NN, dots, 
                   &one, &doubleone, CC[i], &one);
            DGEMV_F77(CHAR_MACRO(T2[0]), &N, &mm, &minusone, UUblock, &NN, dots,
                   &one, &doubleone, UU[i], &one);
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

      dble_tmp = sqrt(AZ_gdot(N, CC[i], CC[i], proc_config));
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

      if ( (iter%print_freq == 0) && proc == 0)
        (void) fprintf(stdout, "%siter: %4d           residual = %e\n",
                       prefix,iter, scaled_r_norm);

      i++;      /* subspace dim. counter dim(K) = i - 1 */

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
                                          proc_config) == AZ_QUIT)) {

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
    }
  }

  if ( (iter%print_freq != 0) && (proc == 0) && (options[AZ_output] != AZ_none)
       && (options[AZ_output] != AZ_warnings))
    (void) fprintf(stdout, "%siter: %4d           residual = %e\n",
		   prefix,iter, scaled_r_norm);

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

} /* AZ_pgmres */

