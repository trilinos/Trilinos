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

#ifdef HAVE_AZTECOO_TEUCHOS
#include "AztecOO_config.h"
#ifdef AZ_ENABLE_TIMEMONITOR
#include "Teuchos_CTimeMonitor.h"
#endif
#endif

extern int az_iterate_id;

void AZ_pgmres (double b[], double x[],double weight[], int options[],
                double params[], int proc_config[], double status[], AZ_MATRIX *Amat,
                AZ_PRECOND *precond, struct AZ_CONVERGE_STRUCT *convergence_info)

     /*******************************************************************************

  This routine uses Saad's restarted Genralized Minimum Residual method to solve
  the nonsymmetric matrix problem Ax = b.

  IMPORTANT NOTE: While the 2-norm of the gmres residual is available, the
  actual residual is not normally computed as part of the gmres algorithm. Thus,
  if the user uses a convergence condition (see AZ_compute_global_scalars())
  that is based on the 2-norm of the residual there is no need to compute the
  residual (i.e. r_avail = AZ_FALSE). However, if another norm of r is
  requested, AZ_compute_global_scalars() sets r_avail = AZ_TRUE and the
  algorithm computes the residual.

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

  register int k;
  int          i, i1, k1, mm, ii;
  int          N, one = 1, iter = 1, r_avail = AZ_FALSE;
  int          precond_flag, print_freq, proc, kspace, first_time = AZ_TRUE;
  double     **v, **hh, *c, *s, *rs, *dots, *tmp, *temp;
  double       *res, init_time = 0.0;
  double       dble_tmp, r_2norm, first_2norm;
  double       rec_residual, scaled_r_norm, true_scaled_r=0.0;
  double       actual_residual = -1.0;
  double       doubleone = 1.0, minusone = -1.0, *dummy = (double *) 0;
  int          *data_org, str_leng;
  char         label[64],suffix[32], prefix[64];


  /* condition number estimation variables */

  double *svbig, *svsml, *vectmp;
  double  big, cc, dble_tmp1, sestpr, small, ss;
  int     ijob;
  double *hhblock, *vblock;
  int kspace_p1,kspace_p2,N_total;
  int aligned_N_total;
  char *T = "T";
  char *T2 = "N";
  int type_orthog, num_orthog_steps;


  /**************************** execution begins ******************************/

  sprintf(suffix," in gmres%d",options[AZ_recursion_level]);/* set string that will be used */
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
  kspace       = options[AZ_kspace];

  num_orthog_steps = 1;
  if (
      options[AZ_orthog]==AZ_classic ||
      options[AZ_orthog]==AZ_double_modified ||
      options[AZ_orthog]==AZ_double_classic) num_orthog_steps = 2;

  type_orthog = 0;
  if (
      options[AZ_orthog]==AZ_modified ||
      options[AZ_orthog]==AZ_single_modified ||
      options[AZ_orthog]==AZ_double_modified) type_orthog = 1;

  /* Initialize some values in convergence info struct */
  convergence_info->print_info = print_freq;
  convergence_info->iteration = 0;
  convergence_info->sol_updated = 0; /* GMRES seldom updates solution */
  convergence_info->epsilon = params[AZ_tol];

  /* allocate memory for required vectors */

  kspace_p2 = kspace + 2;
  kspace_p1 = kspace + 1;
  N_total   = N + data_org[AZ_N_external] + 1;
  /* +1: make sure everybody allocates something */

  /* Note: temp must be aligned on the Intel Paragon so  */
  /* that the quad load inside the assembly dgemv works. */

  sprintf(label,"general%s",suffix);
  temp   = AZ_manage_memory((3*kspace_p2 + 5*kspace_p1 + N_total +
                             (kspace+1)*kspace_p1)
                            *sizeof(double),AZ_ALLOC,
                            AZ_SYS+az_iterate_id,label, &i);

  dots   = &(temp[  N_total]);
  tmp    = &(dots[  kspace_p2]);
  rs     = &(tmp[   kspace_p2]);
  c      = &(rs[    kspace_p2]);
  s      = &(c[     kspace_p1]);
  svbig  = &(s[     kspace_p1]);
  svsml  = &(svbig[ kspace_p1]);
  vectmp = &(svsml[ kspace_p1]);
  hhblock= &(vectmp[kspace_p1]);

  sprintf(label,"ptrs%s",suffix);
  v     = (double **) AZ_manage_memory(2*kspace_p2*sizeof(double *),AZ_ALLOC,
                                       AZ_SYS+az_iterate_id, label, &i);
  hh   = &(v[kspace_p2]);

  aligned_N_total  = N_total;   /* The vectors in 'v' must be aligned on */
  aligned_N_total += N_total%2; /* the Intel Paragon so that the quad    */
                                /* load inside the assembly dgemv works. */

  sprintf(label,"vblock%s",suffix);
  vblock = AZ_manage_memory((kspace+1)*aligned_N_total*sizeof(double),AZ_ALLOC,
                            AZ_SYS+az_iterate_id,label, &i);

  for (k = 0; k < kspace+1; k++) {
    hh[k] = &(hhblock[k*kspace_p1]);
    v[k]  = &(vblock[k*aligned_N_total]);
  }


  AZ_compute_residual(b, x, v[0], proc_config, Amat);

  /*
   * Compute a few global scalars:
   *     1) ||r||                corresponding to options[AZ_conv]
   *     2) scaled ||r||         corresponding to options[AZ_conv]
   *     3) r_2norm = <r,r>      corresponding to options[AZ_conv]
   */


  /* Change to support AztecOO_StatusTest:
     Even though AZ_compute_global_scalars can compute r_2norm for us,
     we want to compute the initial residual prior to calling 
     this function because we need this value passed in to 
     AZ_compute_global_scalars in order to be consistent with subsequent calls
     to the same function later.  
  */
  r_2norm = DDOT_F77(&N, v[0], &one, v[0], &one);
  AZ_gdot_vec(1, &r_2norm, &rec_residual, proc_config);  
  r_2norm = sqrt(r_2norm);
  rec_residual = r_2norm;

  AZ_compute_global_scalars(Amat, x, b, v[0],
                            weight, &rec_residual, &scaled_r_norm, options,
                            data_org, proc_config, &r_avail, dummy,
			    dummy, dummy, convergence_info);
  true_scaled_r = scaled_r_norm;

  r_2norm   = rec_residual;

  if (r_avail) {
    sprintf(label,"res%s",suffix);
    res = AZ_manage_memory(N_total*sizeof(double),AZ_ALLOC,
                           AZ_SYS+az_iterate_id,label,&i);
  }
  else res = (double *) NULL;


  if ( (options[AZ_output] != AZ_none) && 
       (options[AZ_output] != AZ_last) &&
       (options[AZ_output] != AZ_warnings) &&
       (options[AZ_output] != AZ_summary) &&
       (options[AZ_conv]!=AZTECOO_conv_test) && (proc == 0) )
    (void) AZ_printf_out("%siter:    0           residual = %e\n",prefix,scaled_r_norm);

  iter = 0;
  while (!(convergence_info->converged) && iter < options[AZ_max_iter] && !(convergence_info->isnan)) {
    convergence_info->iteration = iter;
    if (r_avail) DCOPY_F77(&N, v[0], &one, res, &one);

    /* v1 = r0/beta */

    dble_tmp    = 1.0 / r_2norm;
    first_2norm = r_2norm;
    DSCAL_F77(&N, &dble_tmp, v[0], &one);

    rs[0] = r_2norm;  /* initialize 1st rhs term of H system */
    i     = 0;

    while (i < kspace && !(convergence_info->converged) && iter < options[AZ_max_iter]
	   && !(convergence_info->isnan)) {
      iter++;
    convergence_info->iteration = iter;
      i1 = i + 1;

      /* v_i+1 = A M^-1 v_i */

      DCOPY_F77(&N, v[i], &one, temp, &one);

      if (iter == 1) init_time = AZ_second();

      if (precond_flag) {

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
        /* Start timer. */
      static int precID = -1;
      precID = Teuchos_startTimer( "AztecOO: Operation Prec*x", precID );
#endif
#endif
        precond->prec_function(temp,options,proc_config,params,Amat,precond);

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Stop timer. */
      Teuchos_stopTimer( precID );
#endif
#endif
      }

      if (iter == 1) status[AZ_first_precond] = AZ_second() - init_time;

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Start timer. */
      static int matvecID = -1;
      matvecID = Teuchos_startTimer( "AztecOO: Operation Op*x", matvecID );
#endif
#endif

      Amat->matvec(temp, v[i1], Amat, proc_config);

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Stop timer. */
      Teuchos_stopTimer( matvecID );
#endif
#endif
      /* Use ||v[i1]|| as estimate for ||A|| in checks below for breakdown. */

      /* Gram-Schmidt orthogonalization */

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Start the timer. */
      static int orthoID = -1;
      orthoID = Teuchos_startTimer( "AztecOO: Orthogonalization", orthoID );
#endif
#endif

      if (type_orthog==0) { /* Classical. Actually, we do */
	                    /* this twice. DGKS method */

        for (k = 0; k <= i; k++) hh[k][i] = 0.0;
        for (ii = 0 ; ii < num_orthog_steps; ii++ ) {
          if (N == 0) for (k = 0; k <= i; k++) dots[k] = 0.0;
          dble_tmp = 0.0; mm = i+1;
#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Start the timer. */
      static int orthoInnerProdID = -1;
      orthoInnerProdID = Teuchos_startTimer( "AztecOO: Ortho (Inner Product)", orthoInnerProdID );
#endif
#endif
          DGEMV_F77(CHAR_MACRO(T[0]), &N, &mm, &doubleone, vblock, &aligned_N_total,
                    v[i1], &one, &dble_tmp, dots, &one);
          AZ_gdot_vec(i1, dots, tmp, proc_config);
#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      Teuchos_stopTimer( orthoInnerProdID );
#endif
#endif
          for (k = 0; k <= i; k++) hh[k][i] += dots[k];

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Start the timer. */
      static int orthoUpdateID = -1;
      orthoUpdateID = Teuchos_startTimer( "AztecOO: Ortho (Update)", orthoUpdateID );
#endif
#endif
          DGEMV_F77(CHAR_MACRO(T2[0]), &N, &mm, &minusone, vblock, &aligned_N_total,
                    dots, &one, &doubleone, v[i1], &one);
#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      Teuchos_stopTimer( orthoUpdateID );
#endif
#endif
        }
      }
      else {                    /* modified */
        for (k = 0; k <= i; k++) hh[k][i] = 0.0;
        for (ii = 0 ; ii < num_orthog_steps; ii++ ) {
	  for (k = 0; k <= i; k++) {
	    hh[k][i] += dble_tmp = AZ_gdot(N, v[k], v[i1], proc_config);
	    dble_tmp = -dble_tmp;
	    DAXPY_F77(&N, &dble_tmp, v[k], &one, v[i1], &one);
	  }
	}
      }

      /* normalize vector */

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      static int orthoNormID = -1;
      orthoNormID = Teuchos_startTimer( "AztecOO: Ortho (Norm)", orthoNormID );
#endif
#endif
      hh[i1][i] = dble_tmp = sqrt(AZ_gdot(N, v[i1], v[i1], proc_config));
#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      Teuchos_stopTimer( orthoNormID );
#endif
#endif
      /* Relative size doesn't matter here, it just needs to be nonzero.
         Some very, very tiny number, such as DBL_MIN, might cause trouble,
         so check for that.
       */
      if (dble_tmp > 100.0 * DBL_MIN )
        dble_tmp  = 1.0 / dble_tmp;
      else
        dble_tmp = 0.0;

      DSCAL_F77(&N, &dble_tmp, v[i1], &one);

#ifdef AZ_ENABLE_TIMEMONITOR
#ifdef HAVE_AZTECOO_TEUCHOS
      /* Stop the timer. */
      Teuchos_stopTimer( orthoID );
#endif
#endif

      /* update factorization of hh by plane rotation */

      for (k = 1; k <= i; k++) {
        k1        = k - 1;
        dble_tmp  = hh[k1][i];
        hh[k1][i] =  c[k1]*dble_tmp + s[k1]*hh[k][i];
        hh[k][i]  = -s[k1]*dble_tmp + c[k1]*hh[k][i];
      }

      /* determine next plane rotation */

      dble_tmp = sqrt(hh[i][i] * hh[i][i] + hh[i1][i] * hh[i1][i]);

      /* Estimate condition number of the GMRES */
      /* least-squares problem using ICE.       */

      if (i == 0) {
        big = dble_tmp;
        small = big;
        svbig[0] = doubleone;
        svsml[0] = doubleone;
      }
      else {
        for (k = 0; k < i; k++) vectmp[k] = hh[k][i];
        vectmp[i] = dble_tmp;
        ijob = 1;
        AZ_DLAIC1_F77(&ijob, &i, svbig, &big, vectmp, &vectmp[i], &sestpr, &ss, &cc);
        big = sestpr;
        DSCAL_F77(&i, &ss, svbig, &one);
        svbig[i] = cc;
        ijob = 2;
        AZ_DLAIC1_F77(&ijob, &i, svsml, &small, vectmp, &vectmp[i], &sestpr, &ss,
                   &cc);
        small = sestpr;
        DSCAL_F77(&i, &ss, svsml, &one);
        svsml[i] = cc;
      }

      /* Check dble_tmp relative to 1.0 here.  Perhaps
         ||A|| would be better.  ||v[i1]|| before it is orthogonalized
         could be used as an estimate of ||A||. */
      if ((small == 0.0) || (dble_tmp  < DBL_EPSILON) ||
          (big/small > params[AZ_ill_cond_thresh]) ) {
        /* (big/small > 1.0e+11) ) {  This is now a parameter */

        /* take most recent solution and get out */

        for (k = 0; k < i1; k++) tmp[k] = rs[k];
        AZ_get_x_incr(options, data_org, proc_config, params, i, hh, tmp,
                      temp, v, Amat, precond, x, &first_time, &(convergence_info->converged), kspace);

        AZ_scale_true_residual(x, b,
                               v[kspace], weight, &actual_residual,
                               &true_scaled_r, options, data_org, proc_config,
                               Amat, convergence_info);

        /* Use same check as above. */
        if (dble_tmp  < DBL_EPSILON) i = AZ_breakdown;
        else                                   i = AZ_ill_cond;

        AZ_terminate_status_print(i, iter, status, rec_residual, params,
                                  true_scaled_r, actual_residual, options,
                                  proc_config);
        return;
      }

      dble_tmp = 1.0 / dble_tmp;
      c[i]     =  hh[i][i]  * dble_tmp;
      s[i]     =  hh[i1][i] * dble_tmp;
      rs[i1]   = -s[i] * rs[i];
      rs[i]   *=  c[i];

      if (r_avail) {
        dble_tmp  = c[i] * rs[i1];
        dble_tmp1 = s[i] * s[i];

        for (k = 0; k < N; k++)
          res[k] = dble_tmp*v[i1][k] + dble_tmp1*res[k];
      }

      /* determine residual norm & test convergence */

      hh[i][i]     = c[i] * hh[i][i] + s[i] * hh[i1][i];
      r_2norm      = fabs(rs[i1]);
      rec_residual = r_2norm;

      /*
       * Compute a few global scalars:
       *     1) ||r||                corresponding to options[AZ_conv]
       *     2) scaled ||r||         corresponding to options[AZ_conv]
       * NOTE: if r_avail = AZ_TRUE or AZ_FIRST is passed in, we perform
       * step 1), otherwise ||r|| is taken as rec_residual.
       */

      AZ_compute_global_scalars(Amat, x, b,
                                res, weight, &rec_residual, &scaled_r_norm,
                                options, data_org, proc_config, &r_avail, dummy,
                                dummy, dummy, convergence_info);

      if ( (iter%print_freq == 0) &&
           (options[AZ_conv]!=AZTECOO_conv_test) && proc == 0)
        (void) AZ_printf_out("%siter: %4d           residual = %e\n",prefix,iter,
                       scaled_r_norm);

      i++;      /* subspace dim. counter dim(K) = i - 1 */

      if (convergence_info->isnan) {
        AZ_terminate_status_print(AZ_breakdown, iter, status, rec_residual, params,
                                  true_scaled_r, actual_residual, options,
                                  proc_config);
        return;
      }
      if ( (i == kspace) || convergence_info->converged || iter == options[AZ_max_iter]) {

        /* update x and set temp to delta x */

        for (k = 0; k <= i1; k++) tmp[k] = rs[k];

        AZ_get_x_incr(options, data_org, proc_config, params, i, hh, tmp,
                      temp, v, Amat, precond,  x, &first_time, &(convergence_info->converged), kspace);

      }

      if (convergence_info->converged) {

        /* compute true residual using 'v[kspace]' as a temporary vector */

        AZ_scale_true_residual(x, b,
                               v[kspace], weight, &actual_residual,
                               &true_scaled_r, options, data_org, proc_config,
                               Amat, convergence_info);


        if (!(convergence_info->converged) && options[AZ_conv]!=AZTECOO_conv_test) {

	  if (AZ_get_new_eps(&(convergence_info->epsilon), scaled_r_norm,
                                          true_scaled_r,
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

        /* restore previous solution */

        if ( (!(convergence_info->converged)) && (i != kspace) )
          for (k = 0; k < N; k++) x[k] = x[k] - temp[k];
      }

      if ( (i == kspace) && !(convergence_info->converged)) {
        if (r_avail)
          for (k = 0; k < N; k++) v[0][k] = res[k];
        else {
          Amat->matvec(temp, v[kspace], Amat, proc_config);

          DSCAL_F77(&N, &first_2norm, v[0], &one);
          DAXPY_F77(&N, &minusone, v[kspace], &one, v[0], &one);
        }
      }
    }
  }

  if ( (iter%print_freq != 0) && (proc == 0) && (options[AZ_output] != AZ_none)
       && (options[AZ_output] != AZ_warnings) && (options[AZ_conv]!=AZTECOO_conv_test))
    (void) AZ_printf_out("%siter: %4d           residual = %e\n",prefix,iter,
                   scaled_r_norm);


  /* check if we exceeded maximum number of iterations */

  if (convergence_info->converged) {
    i = AZ_normal;
    scaled_r_norm = true_scaled_r;
  }
  else if (convergence_info->isnan) i = AZ_breakdown;
  else i = AZ_maxits;

  AZ_terminate_status_print(i, iter, status, rec_residual, params,
                            scaled_r_norm, actual_residual, options,
                            proc_config);

} /* AZ_pgmres */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_get_x_incr(int options[], int data_org[], int
                   proc_config[], double params[], int i, double **hh, double
                   *rs, double *temp, double **v, AZ_MATRIX *Amat, AZ_PRECOND
                   *precond, double x[], int *first_time, int *converged,
                   int kspace)

     /*******************************************************************************

  This routine is normally invoked from GMRES and is used to compute the
  increment to the solution (as well as the new solution) that should be applied as a
  result of solving the upper hessenberg system produces by GMRES.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               On input, contains the initial guess. On output contains the
                   solution to the linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  params:          Drop tolerance and convergence tolerance info.

  i:               The size of the Krylov subspace that has be generated.

  hh:              The upper triangular matrix produced by gmres (after reducing
                   the upper hessenberg matrix) after sweeping out i steps.

  rs:              The projected right hand side produced by gmres.

  temp:            vector of length data_org[AZ_num_int_unk] +
                   data_org[AZ_num_bord_unk].  On output to AZ_get_x_incr(),
                   temp contains the increment that must be added the current
                   approximation to get the new gmres approximate solution.

  v:               Orthogonal vectors produced by Gram-Schmidt process that
                   span the Krylov space swept out.


  Amat:            Structure used to represent the matrix (see az_aztec.h
                   and Aztec User's Guide).

  precond:         Structure used to represent the preconditionner
                   (see az_aztec.h ad Aztec User's Guide).


*******************************************************************************/

{

  /* local variables */

  int precond_flag;
  int    ii, k, k1, j, N;
  double t, doubleone = 1.0, update_norm = 1.0;
  int    one = 1;

  /**************************** execution begins ******************************/

  if (i == 0) return;

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  /* solve upper triangular system, compute solution */

  i--;  /* set i = Krylov subspace dimension */

  rs[i] /= hh[i][i];
  for (ii = 1; ii <= i; ii++) {
    k  = i - ii;
    k1 = k + 1;
    t  = rs[k];
    for (j = k1; j <= i; j++) t -= hh[k][j] * rs[j];
    rs[k] = t / hh[k][k];
  }

  /*
   * done with back substitution; form linear combination to get solution
   */

  precond_flag = options[AZ_precond];
  if (options[AZ_check_update_size] & *converged) {
    for (j = 0; j < N; j++) temp[j] = v[i][j];
    if (precond_flag) precond->prec_function(temp,options,proc_config,params, Amat,precond);
    update_norm = fabs(rs[i]*sqrt(AZ_gdot(N, temp, temp, proc_config)));
  }

  for (j = 0; j < N; j++) temp[j] = 0.0;

  for (j = 0; j <= i; j++) {
    t = rs[j];
    DAXPY_F77(&N, &t, v[j], &one, temp, &one);
  }

  if (precond_flag) precond->prec_function(temp,options,proc_config,params,
                                           Amat,precond);

  DAXPY_F77(&N, &doubleone, temp, &one, x, &one);

  if (options[AZ_check_update_size] & *converged) {
    *converged =AZ_compare_update_vs_soln(N,update_norm, rs[i],temp,x,
                                          params[AZ_update_reduction],options[AZ_output],proc_config,first_time);

    /* restore previous solution */

    if ( (!(*converged)) && (i != kspace) ) {
      doubleone = -1.;
      DAXPY_F77(&N, &doubleone, temp, &one, x, &one);
    }
  }

} /* AZ_get_x_incr */
