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
#include "az_lapack_wrappers.h"

extern int az_iterate_id;

static void compute_condnum_tridiag_sym( int N, double *diag, double *offdiag,
					 char prefix[], int options[],
					 int proc_config[],
                                         double*, double*, 
					 double *ConditionNumber );

void AZ_pcg_f_condnum(double b[], double x[], double weight[], int options[],
			     double params[], int proc_config[],double status[],
			     AZ_MATRIX *Amat, AZ_PRECOND *precond,
			     struct AZ_CONVERGE_STRUCT *convergence_info )

     /*******************************************************************************

  Conjugate Gradient algorithm to solve the symmetric matrix problem Ax = b.

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

  precond:         Structure used to represent the preconditioner
                   (see file az_aztec.h and Aztec User's Guide).
*******************************************************************************/



{

  /* local variables */

  register int i;
  int          N, NN, one = 1, iter = 1, r_avail = AZ_TRUE, j;
  int          precond_flag, print_freq, proc, brkdown_will_occur = AZ_FALSE;
  double       alpha, beta = 0.0, nalpha, true_scaled_r=-1.0;
  double      *r, *z, *p, *ap, actual_residual = -1.0;
  double       r_z_dot, r_z_dot_old, p_ap_dot, rec_residual=-1.0;
  double       scaled_r_norm=-1.0, brkdown_tol = DBL_EPSILON;
  int          *data_org, str_leng, first_time = AZ_TRUE;
  char         label[64],suffix[32], prefix[64];

  double **saveme, *ptap;
  int *kvec_sizes = NULL, current_kept = 0;
  double *dots;
  double doubleone = 1., dzero = 0.;
  char *T = "T";
  char *T2 = "N";
  double *block;

  /* condition number estimate from Lanczos matrix */
  double beta_old = 0, p_ap_dot_old = 0;
  int N_lanczos_max = options[AZ_max_iter];
  double * diag_T = NULL;
  double * offdiag_T = NULL;
  int N_lanczos=0;
  double r_z_dot_old2 = 0;
  double lambda_min, lambda_max;
  double ConditionNumber;

  /**************************** execution begins ******************************/

  diag_T = (double *) malloc( sizeof(double) * N_lanczos_max ); 
  offdiag_T = (double *) malloc( sizeof(double) * (N_lanczos_max-1) ); 

  sprintf(suffix," in cg%d",options[AZ_recursion_level]);  /* set string that will be used */
                                                           /* for manage_memory label      */
  /* set prefix for printing */

  str_leng = 0;
  for (i = 0; i < 16; i++) prefix[str_leng++] = ' ';
  for (i = 0 ; i < options[AZ_recursion_level]; i++ ) {
    prefix[str_leng++] = ' '; prefix[str_leng++] = ' '; prefix[str_leng++] = ' ';
    prefix[str_leng++] = ' '; prefix[str_leng++] = ' ';
  }
  prefix[str_leng] = '\0';


  /* pull needed values out of parameter arrays */

  data_org = Amat->data_org;

  N            = data_org[AZ_N_internal] + data_org[AZ_N_border];


  precond_flag = options[AZ_precond];
  proc         = proc_config[AZ_node];
  print_freq   = options[AZ_print_freq];

  /* Initialize some values in convergence info struct */
  convergence_info->print_info = print_freq;
  convergence_info->iteration = 0;
  convergence_info->sol_updated = 1; /* CG always updates solution */
  convergence_info->epsilon = params[AZ_tol]; /* Test against this */

  /* allocate space for necessary vectors */

  NN = N + data_org[AZ_N_external];
  if (NN == 0) NN++;  /* make sure everybody allocates something */
  NN = NN + (NN%2);   /* make sure things are aligned for assembly */
                      /* matvec on paragon. */



  sprintf(label,"z%s",suffix);
  p  = (double *) AZ_manage_memory(4*NN*sizeof(double),AZ_ALLOC,
                                   AZ_SYS+az_iterate_id, label, &j);
  r  = &(p[1*NN]);
  z  = &(p[2*NN]);
  ap = &(p[3*NN]);

  AZ_compute_residual(b, x, r, proc_config, Amat);

  if (options[AZ_apply_kvecs]) {
    AZ_compute_global_scalars(Amat, x, b, r,
                              weight, &rec_residual, &scaled_r_norm, options,
                              data_org, proc_config, &r_avail,NULL, NULL, &r_z_dot,
                              convergence_info);
    AZ_space_for_kvecs(AZ_OLD_ADDRESS, &kvec_sizes, &saveme,
                       &ptap, options, data_org, suffix,
                       proc_config[AZ_node], &block);
    dots = (double *) AZ_allocate(2*kvec_sizes[AZ_Nkept]*sizeof(double));
    if (dots == NULL) {
      printf("Not space to apply vectors in CG\n");
      exit(1);
    }
    DGEMV_F77(CHAR_MACRO(T[0]),&N,&(kvec_sizes[AZ_Nkept]),&doubleone,block,&N, r, &one, &dzero, dots, &one);
    AZ_gdot_vec(kvec_sizes[AZ_Nkept], dots, &(dots[kvec_sizes[AZ_Nkept]]), proc_config);
    for (i = 0; i < kvec_sizes[AZ_Nkept]; i++) dots[i] = dots[i]/ptap[i];
    DGEMV_F77(CHAR_MACRO(T2[0]), &N, &(kvec_sizes[AZ_Nkept]), &doubleone, block, &N, dots, &one, &doubleone,
              x,  &one);

    AZ_free(dots);
    AZ_compute_residual(b, x, r, proc_config, Amat);
    if ((options[AZ_output] != AZ_none) && (proc == 0))
      printf("\t\tApplied Previous Krylov Vectors ... \n\n");
  }
  if (options[AZ_keep_kvecs] > 0)
    AZ_space_for_kvecs(AZ_NEW_ADDRESS, &kvec_sizes, &saveme,
                       &ptap, options, data_org, suffix,
                       proc_config[AZ_node], &block);



  /*  z = M r */
  /*  p = 0   */

  DCOPY_F77(&N, r, &one, z, &one);
  status[AZ_first_precond] = AZ_second();
  if (precond_flag)
    precond->prec_function(z,options,proc_config,params,Amat,precond);

  status[AZ_first_precond] = AZ_second() - status[AZ_first_precond];

  for (i = 0; i < N; i++ ) p[i] = 0.0;

  /* compute a few global scalars:                                 */
  /*     1) ||r||                corresponding to options[AZ_conv] */
  /*     2) scaled ||r||         corresponding to options[AZ_conv] */
  /*     3) r_z_dot = <z, r>                                       */

  AZ_compute_global_scalars(Amat, x, b, r,
                            weight, &rec_residual, &scaled_r_norm, options,
                            data_org, proc_config, &r_avail,r, z, &r_z_dot,
                            convergence_info);
  true_scaled_r = scaled_r_norm;

  if ((options[AZ_output] != AZ_none) &&
      (options[AZ_output] != AZ_last) &&
      (options[AZ_output] != AZ_warnings) &&
      (options[AZ_output] != AZ_summary) &&
      (options[AZ_conv]!=AZTECOO_conv_test) && (proc == 0))
    {
      (void) AZ_printf_out("%siter:    0           residual = %e\n",
                     prefix,scaled_r_norm);
      AZ_flush_out();
    }


  for (iter = 1; iter <= options[AZ_max_iter] && !(convergence_info->converged) && 
	 !(convergence_info->isnan); iter++ ) {
    convergence_info->iteration = iter;

    /* p  = z + beta * p */
    /* ap = A p          */

    for (i = 0; i < N; i++) p[i] = z[i] + beta * p[i];
    Amat->matvec(p, ap, Amat, proc_config);

    if ((options[AZ_orth_kvecs]) && (kvec_sizes != NULL)) {
      for (i = 0; i < current_kept; i++) {
        alpha = -AZ_gdot(N, ap, saveme[i], proc_config)/ptap[i];
        DAXPY_F77(&N, &alpha,  saveme[i],  &one, p, &one);
      }
      if (current_kept > 0) Amat->matvec(p, ap, Amat, proc_config);
    }

    p_ap_dot = AZ_gdot(N, p, ap, proc_config);
    if (fabs(p_ap_dot) < brkdown_tol) {

      /* possible problem */

      if (AZ_breakdown_f(N, p, ap, p_ap_dot, proc_config)) {

        /* something wrong */

        AZ_scale_true_residual(x, b, ap,
                               weight, &actual_residual, &true_scaled_r,
                               options, data_org, proc_config, Amat,
                               convergence_info);
        AZ_terminate_status_print(AZ_breakdown, iter, status, rec_residual,
                                  params, true_scaled_r, actual_residual,
                                  options, proc_config);

	compute_condnum_tridiag_sym( N_lanczos-2, diag_T, offdiag_T, prefix,
				     options, proc_config, &lambda_min,
                                     &lambda_max, &ConditionNumber);
	status[AZ_lambda_min] = lambda_min;
	status[AZ_lambda_max] = lambda_max;
	status[AZ_condnum] = ConditionNumber;
        free((void*)diag_T);     diag_T = NULL;
        free((void*)offdiag_T);  offdiag_T = NULL;
        return;
      }
      else brkdown_tol = 0.1 * fabs(p_ap_dot);
    }

    alpha  = r_z_dot / p_ap_dot;
    nalpha = -alpha;

    /* x = x + alpha*p  */
    /* r = r - alpha*Ap */
    /* z = M^-1 r       */

    DAXPY_F77(&N, &alpha,  p,  &one, x, &one);

    if (iter <= options[AZ_keep_kvecs]) {
      DCOPY_F77(&N, p, &one, saveme[iter-1], &one);
      ptap[iter-1] = p_ap_dot ;
      kvec_sizes[AZ_Nkept]++;
      current_kept = kvec_sizes[AZ_Nkept];
    }
    /*
      else {
      i = (iter-1)%options[AZ_keep_kvecs];
      DCOPY_F77(&N, p, &one, saveme[i], &one);
      ptap[i] = p_ap_dot ;
      }
    */
    DAXPY_F77(&N, &nalpha, ap, &one, r, &one);
    DCOPY_F77(&N, r, &one, z, &one);

    if (precond_flag) precond->prec_function(z,options,proc_config,params,Amat,precond);

    r_z_dot_old = r_z_dot;

    /* compute a few global scalars:                                 */
    /*     1) ||r||                corresponding to options[AZ_conv] */
    /*     2) scaled ||r||         corresponding to options[AZ_conv] */
    /*     3) r_z_dot = <z, r>                                       */

    AZ_compute_global_scalars(Amat, x, b, r,
                              weight, &rec_residual, &scaled_r_norm, options,
                              data_org, proc_config, &r_avail, r, z, &r_z_dot,
                              convergence_info);

    if (brkdown_will_occur) {
      AZ_scale_true_residual( x, b, ap,
                              weight, &actual_residual, &true_scaled_r, options,
                              data_org, proc_config, Amat,convergence_info);
      AZ_terminate_status_print(AZ_breakdown, iter, status, rec_residual,
                                params, true_scaled_r, actual_residual, options,
                                proc_config);

      compute_condnum_tridiag_sym( N_lanczos-2, diag_T, offdiag_T, prefix,
				   options, proc_config, &lambda_min,
                                   &lambda_max, &ConditionNumber);
      status[AZ_lambda_min] = lambda_min;
      status[AZ_lambda_max] = lambda_max;
      status[AZ_condnum] = ConditionNumber;
      free((void*)diag_T);    diag_T = NULL;
      free((void*)offdiag_T); offdiag_T = NULL;

      return;
    }

    beta = r_z_dot / r_z_dot_old;

    if (fabs(r_z_dot) < brkdown_tol) {

      /* possible problem */

      if (AZ_breakdown_f(N, r, z, r_z_dot, proc_config))
        brkdown_will_occur = AZ_TRUE;
      else
        brkdown_tol = 0.1 * fabs(r_z_dot);
    }

    if ( (iter%print_freq == 0) && (options[AZ_conv]!=AZTECOO_conv_test) && proc == 0 )
      {
        (void) AZ_printf_out("%siter: %4d           residual = %e\n", prefix, iter,
                       scaled_r_norm);
        AZ_flush_out();
      }

    /* convergence tests */

    if (options[AZ_check_update_size] & convergence_info->converged)
      convergence_info->converged = AZ_compare_update_vs_soln(N, -1.,alpha, p, x,
							      params[AZ_update_reduction],
							      options[AZ_output], proc_config, &first_time);


    if (convergence_info->converged) {
      AZ_scale_true_residual(x, b, ap,
                             weight, &actual_residual, &true_scaled_r, options,
                             data_org, proc_config, Amat, convergence_info);
      
      
      
      /*
       * Note: epsilon and params[AZ_tol] may not be equal due to a previous
       * call to AZ_get_new_eps().
       */
      
      if (!(convergence_info->converged) && options[AZ_conv]!=AZTECOO_conv_test) {

	if (AZ_get_new_eps(&(convergence_info->epsilon), scaled_r_norm, true_scaled_r,
			   options, proc_config) == AZ_QUIT) {

	  /*
	   * Computed residual has converged, actual residual has not converged,
	   * AZ_get_new_eps() has decided that it is time to quit.
	   */
	  
	  AZ_terminate_status_print(AZ_loss, iter, status, rec_residual, params,
				    true_scaled_r, actual_residual, options,
				    proc_config);

	  /* compute the eigenvalues of the Lanczos matrix */
	  compute_condnum_tridiag_sym( N_lanczos-2, diag_T, offdiag_T, prefix,
				       options, proc_config, &lambda_min,
                                       &lambda_max, &ConditionNumber);
	  status[AZ_lambda_min] = lambda_min;
	  status[AZ_lambda_max] = lambda_max;
	  status[AZ_condnum] = ConditionNumber;
	  free((void*)diag_T);     diag_T = NULL;
	  free((void*)offdiag_T);  offdiag_T = NULL;
	
	  return;
	}
      }
    }

    /* estimation of the condition number.
       Fill the tridiagonal Lanczos matrix, stored in the arrays
       diag_T and offdiag_T */

    if( iter>1 ) {

      diag_T[N_lanczos]    = (pow(beta_old,2) * p_ap_dot_old + p_ap_dot) / r_z_dot_old;
      offdiag_T[N_lanczos-1] = -beta_old * p_ap_dot_old / (sqrt( r_z_dot_old * r_z_dot_old2));

      N_lanczos++;
      
    } else {

      diag_T[N_lanczos]    = p_ap_dot / r_z_dot_old;
      N_lanczos++;
    }

    r_z_dot_old2 = r_z_dot_old;
      
    beta_old = beta;
    p_ap_dot_old = p_ap_dot;

    /* up to here */
    
  }
  iter--;
  if ( (iter%print_freq != 0) && (proc == 0) && (options[AZ_output] != AZ_none)
       && (options[AZ_output] != AZ_warnings) &&
      (options[AZ_conv]!=AZTECOO_conv_test) )
    {
      (void) AZ_printf_out("%siter: %4d           residual = %e\n", prefix, iter,
                     scaled_r_norm);
      AZ_flush_out();
    }

  /* check if we exceeded maximum number of iterations */

  if (convergence_info->converged) {
    i = AZ_normal; scaled_r_norm = true_scaled_r; }
  else if (convergence_info->isnan) i = AZ_breakdown;
  else
    i = AZ_maxits;

  AZ_terminate_status_print(i, iter, status, rec_residual, params,
                            scaled_r_norm, actual_residual, options,
                            proc_config);

  /* compute the condition number */
  
  compute_condnum_tridiag_sym( N_lanczos-2, diag_T, offdiag_T, prefix,
			       options, proc_config, &lambda_min,
                               &lambda_max, &ConditionNumber);
  status[AZ_lambda_min] = lambda_min;
  status[AZ_lambda_max] = lambda_max;
  status[AZ_condnum] = ConditionNumber;
	  
  if( diag_T != NULL ) {
    free((void *)diag_T); diag_T = NULL;
  }
  if( offdiag_T != NULL ) {
    free((void *)offdiag_T); offdiag_T = NULL;
  }
} /* AZ_pcg_condnum */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void compute_condnum_tridiag_sym( int N, double *diag, double *offdiag,
					 char prefix[], int options[],
					 int proc_config[],
                                         double* lambda_min,
                                         double* lambda_max,
					 double *ConditionNumber )

  /* diag ==      double vector of size N, containing the diagonal
                  elements of A
     offdiag ==   double vector of size N-1, containing the offdiagonal
                  elements of A. Note that A is supposed to be symmatric
  */		  
{
  
  int N_split;
  double zero_double = 0.0;

  double * eigenvalues, * work;
  int * iwork, info, * iblock, * isplit;
  int dummy, N_eigs;
  double double_dummy;
  char char_A = 'A';
  char char_E = 'E';
  double smallest, largest;

  eigenvalues = (double *) malloc( sizeof(double) * N);
  work        = (double *) malloc( sizeof(double) * N * 4 );
  iwork       = (int *) malloc( sizeof(int) * N * 3 );
  iblock      = (int *) malloc( sizeof(int) * N );
  isplit      = (int *) malloc( sizeof(int) * N );

  if( N > 2 ) {
    
    DSTEBZ_F77( CHAR_MACRO(char_A), CHAR_MACRO(char_E),
		&N, &double_dummy, &double_dummy,
		&dummy, &dummy,
		&zero_double, diag, offdiag, &N_eigs, &N_split,
		eigenvalues, iblock, isplit, work, iwork, &info );

    smallest = eigenvalues[0];
    largest = eigenvalues[N_eigs-1];

  } else {
    if( proc_config[AZ_node] == 0 )
      printf( "\n%sWarning : The Lanczos matrix is too small\n",
              prefix );
    
    smallest = 1.0;
    largest = 1.0;
    
  }
  
  /* Simply a check that all the eigenvalues are positive */
  if( smallest <= 0.0 ) {
    if( proc_config[AZ_node] == 0 )
      printf( "\n%sWarning : The smallest eigenvalue of the Lanczos matrix\n"
              "%sis negative or zero (%e)\n",
              prefix,
              prefix,
              smallest );
    
  } 
    
  *lambda_min = smallest;
  *lambda_max = largest;
  *ConditionNumber = largest/smallest;
  
  free((void *) eigenvalues );
  free((void *) work );
  free((void *) iwork );
  free((void *) iblock );
  free((void *) isplit );

  if( proc_config[AZ_node] != 0 ) return;
  
  if ( (options[AZ_output] != AZ_none) && (options[AZ_output] != AZ_last) && (options[AZ_output] != AZ_summary) &&
       (options[AZ_output] != AZ_warnings) && (proc_config[AZ_node] == 0) ) {

    printf("\n%s-----------------------------------------------------\n",
	   prefix);
    
    printf("\n%sAnalysis of the Lanczos matrix of\n"
	   "%sthe preconditioned system:\n\n"
	   "%ssmallest eigenvalue          = %e\n"
	   "%slargest eigenvalue           = %e\n"
	   "\n%sestimated condition number   = %e\n",
	   prefix, 
	   prefix, prefix, smallest, prefix, largest,
	   prefix,  *ConditionNumber );

    printf("\n%s-----------------------------------------------------\n",
	   prefix);
  }  
  
  return;
  
} /* compute_condnum_tridiag_sym */
