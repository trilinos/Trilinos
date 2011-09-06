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
/*---------------- External Definitions -------------------------------------*/

#ifdef eigen
extern void AZ_do_Jacobi(double val[], int indx[], int bindx[], int rpntr[],
                     int cpntr[], int bpntr[], double x[], double b[],
                     double temp[], int options[], int data_org[],
                     int proc_config[], double params[], int flag);
#endif
extern void AZ_calc_blk_diag_inv(double *val, int *indx, int *bindx, int *rpntr,
                     int *cpntr, int *bpntr, double *d_inv, int *d_indx,
                     int *d_bindx, int *d_rpntr, int *d_bpntr,
                     int data_org[]);
extern void jacobi(double val[], double x[], int data_org[]);

extern int az_iterate_id;

/*---------------------------------------------------------------------------*/


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_precondition(double x[], int input_options[], int proc_config[],
                     double input_params[], AZ_MATRIX *Amat, 
		     AZ_PRECOND *input_precond)


/*******************************************************************************

  This routine calls appropriate sparse matrix preconditioner.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============


  x:               On input, contains the current solution. On output contains
                   the preconditioned solution to the linear system.

  options:         Determines specific solution method and other parameters.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  params:          Drop tolerance and convergence tolerance info.

  Amat:            Structure used to represent the matrix (see az_aztec.h
                   and Aztec User's Guide).

  precond:         Structure used to represent the preconditionner
                   (see file az_aztec.h and Aztec User's Guide).

 * --------------------------------------------------------------------

 Related routines:

   scaling routines:
        AZ_block_diagonal_scaling -- block-diagonally scales sparse matrix
                                     problem.
        AZ_row_sum_scaling        -- row sum scales sparse matrix problem.
        sym_diagonal_scaling      -- diagonaly scales symm. sparse problem.
        sym_row_sum_scaling       -- row sum scales symmetric sparse problem.

   preconditioners:
        jacobi                 -- point Jacobi method.
        AZ_polynomial_expansion-- Polynomial expansion; Neumann series and
                                  least squares.
        domain decomposition   -- Block solvers (LU , ILU or ILUT) used on 
                                  each processor. The blocks are either
                                  non-overlapping or overlapping.
        icc                    -- incomplete sparse Choleski (symmetric
                                  version).

*******************************************************************************/

{

  /* local variables */

  int            ione = 1;
  double        *temp;
  int            m, N, k, length;
  int            i, step, j;
  static int    *d2_indx,*d2_bindx,*d2_rpntr,*d2_bpntr;
  static double *d2_inv;
  static AZ_MATRIX *Dmat;
  int            tsize, multilevel_flag = 0, max_externals;
  static int     previous_factors = -1;
  double        *v, *y;
  char          *yo = "precond: ";
  int          *data_org, *bindx, *indx, *cpntr, *rpntr, *bpntr;
  double       *val;
  char         label[64],suffix[32];
  char         tag[80];
  double       *current_rhs, *orig_rhs = NULL, *x_precond = NULL;
  int          *options, *ioptions, N_fixed, *fixed_pts;
  double       *params,  *iparams, *istatus;
  AZ_MATRIX    *Aptr, *Pmat;
  AZ_PRECOND   *Pptr, *precond;
  struct AZ_SCALING *Sptr;
  int          opt_save1, opt_save2, opt_save3, opt_save4, opt_save5, *itemp;
  double       *tttemp, norm1, *dtemp;
#ifdef TIMING
  double       ttt;
#endif


#ifdef eigen
  double         *tb, *tr;
#endif

  /**************************** execution begins ******************************/
#ifdef TIMING
  ttt = AZ_second();
#endif

  precond = input_precond;

  sprintf(suffix," in precond%d",input_options[AZ_recursion_level]);  
                                              /* set string that will be used */
                                              /* for manage_memory label      */

  data_org = precond->Pmat->data_org;
  options  = input_options;
  params   = input_params;

  m    = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];
  N    = data_org[AZ_N_internal] + data_org[AZ_N_border];
  max_externals = Amat->data_org[AZ_N_external];
  if (max_externals < data_org[AZ_N_external]) 
     max_externals = data_org[AZ_N_external];

  current_rhs = x; 
  if (options[AZ_precond] == AZ_multilevel) {

     /* make extra vectors to hold rhs and residual */

     sprintf(tag,"orig_rhs %s",precond->context->tag);
     orig_rhs = AZ_manage_memory((N+max_externals)*sizeof(double),
                               AZ_ALLOC, AZ_SYS+az_iterate_id,tag,&i);
     sprintf(tag,"x_prec %s",precond->context->tag);
     x_precond    = AZ_manage_memory((N+max_externals)*sizeof(double),
                               AZ_ALLOC, AZ_SYS+az_iterate_id, tag,&i);
     for (i = 0 ; i < N; i++) x_precond[i] = 0.0;
     for (i = 0 ; i < N; i++) orig_rhs[i] = current_rhs[i];
     multilevel_flag = 1;
     options = precond->options;
     params  = precond->params;
  }

  do {
     data_org = precond->Pmat->data_org;
     val      = precond->Pmat->val;
     bindx    = precond->Pmat->bindx;
     cpntr    = precond->Pmat->cpntr;
     indx     = precond->Pmat->indx;
     rpntr    = precond->Pmat->rpntr;
     bpntr    = precond->Pmat->bpntr;
     if (max_externals < data_org[AZ_N_external]) 
        max_externals = data_org[AZ_N_external];

     switch (options[AZ_precond]) {
     case AZ_none:
     break;

     case AZ_Jacobi:
        if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
           for (i = 0; i < N; i++) current_rhs[i] /= val[i];

           if (options[AZ_poly_ord] > 1) {
              sprintf(tag,"v_prec %s",precond->context->tag);
              v = AZ_manage_memory((N+max_externals)*sizeof(double),
                                    AZ_ALLOC, AZ_SYS+az_iterate_id, tag, &i);
              sprintf(tag,"y_prec %s",precond->context->tag);
              y = AZ_manage_memory(N*sizeof(double), AZ_ALLOC, AZ_SYS+az_iterate_id, tag,&i);
              for (i = 0; i < N; i++) v[i] = current_rhs[i];

              for (step = 1; step < options[AZ_poly_ord]; step++) {
                 Amat->matvec(v, y, Amat, proc_config);
                 for(i = 0; i < N; i++) v[i] += current_rhs[i] - y[i] / val[i];
              }
              for (i = 0; i < N; i++) current_rhs[i] = v[i];
           }
        }
        else if (data_org[AZ_matrix_type] == AZ_USER_MATRIX) {
           if (options[AZ_pre_calc] < AZ_sys_reuse) {
              sprintf(tag,"d2_inv %s",precond->context->tag);
              d2_inv   = (double *) AZ_manage_memory(N*sizeof(double),AZ_ALLOC,
						data_org[AZ_name],tag,&i);
              Pmat = precond->Pmat;
              if ( (Pmat->N_nz < 0) || (Pmat->max_per_row < 0)) 
                 AZ_matfree_Nnzs(Pmat);

              if ( (Pmat->getrow == NULL) && (N != 0) ) {
                 AZ_printf_err("Error: Only matrices with getrow() defined via ");
                 AZ_printf_err("AZ_set_MATFREE_getrow(...) can do Jacobi preconditioning\n");
                 exit(1);
              }
              sprintf(tag,"dtemp %s",precond->context->tag);
              dtemp = (double *) AZ_manage_memory(Pmat->max_per_row*
				                sizeof(double),AZ_ALLOC,
						data_org[AZ_name],tag,&i);
              sprintf(tag,"itemp %s",precond->context->tag);
              itemp = (int    *) AZ_manage_memory(Pmat->max_per_row*
				                sizeof(int   ),AZ_ALLOC,
						data_org[AZ_name],tag,&i);
  
	      for (i = 0; i < N; i++) {
                 Pmat->getrow(itemp,dtemp,&length,Pmat,1,&i,Pmat->max_per_row);
                 for (k =0; k < length; k++) 
                    if (itemp[k] == i) break;

                 if (k == length) d2_inv[i] = 0.0; /* no diagonal */
                 else d2_inv[i] = 1./dtemp[k];
              }
           }
           for (i = 0; i < N; i++) current_rhs[i] *= d2_inv[i];

           if (options[AZ_poly_ord] > 1) {
              sprintf(tag,"v_prec %s",precond->context->tag);
              v = AZ_manage_memory((N+max_externals)*sizeof(double),
                                    AZ_ALLOC, AZ_SYS+az_iterate_id, tag, &i);
              sprintf(tag,"y_prec %s",precond->context->tag);
              y = AZ_manage_memory(N*sizeof(double), AZ_ALLOC, AZ_SYS+az_iterate_id, tag,&i);
              for (i = 0; i < N; i++) v[i] = current_rhs[i];

              for (step = 1; step < options[AZ_poly_ord]; step++) {
                 Amat->matvec(v, y, Amat, proc_config);
                 for(i = 0; i < N; i++) v[i] += current_rhs[i] - y[i]*d2_inv[i];
              }
              for (i = 0; i < N; i++) current_rhs[i] = v[i];
           }
        }
        else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
           /* block Jacobi preconditioning */

           if (options[AZ_pre_calc] < AZ_sys_reuse) {
              /* First, compute block-diagonal inverse */
              /* (only if not already computed)        */

              tsize = 0;
              for (i = 0; i < m; i++)
                 tsize += (rpntr[i+1] - rpntr[i]) * (cpntr[i+1] - cpntr[i]);

                 sprintf(tag,"d2_indx %s",precond->context->tag);
                 d2_indx  = (int *) AZ_manage_memory((m+1)*sizeof(int),AZ_ALLOC,
                                            data_org[AZ_name], tag, &i);
                 sprintf(tag,"d2_bindx %s",precond->context->tag);
                 d2_bindx = (int *) AZ_manage_memory(m*sizeof(int), AZ_ALLOC,
                                            data_org[AZ_name], tag, &i);
                 sprintf(tag,"d2_rpntr %s",precond->context->tag);
                 d2_rpntr = (int *) AZ_manage_memory((m+1)*sizeof(int),AZ_ALLOC,
                                            data_org[AZ_name], tag, &i);
                 sprintf(tag,"d2_bpntr %s",precond->context->tag);
                 d2_bpntr = (int *) AZ_manage_memory((m+1)*sizeof(int),AZ_ALLOC,
                                            data_org[AZ_name], tag, &i);
                 sprintf(tag,"d2_inv %s",precond->context->tag);
                 d2_inv   = (double *) AZ_manage_memory(tsize*sizeof(double),
                                            AZ_ALLOC, data_org[AZ_name],tag,&i);
                 d2_bpntr[0] = 0;
                 sprintf(tag,"dmat_calk_binv %s",precond->context->tag);
                 Dmat     = (AZ_MATRIX *) AZ_manage_memory(sizeof(AZ_MATRIX), 
                                            AZ_ALLOC,data_org[AZ_name],tag,&i);

                 Dmat->rpntr         = d2_rpntr;   Dmat->cpntr    = d2_rpntr;
                 Dmat->bpntr         = d2_bpntr;   Dmat->bindx    = d2_bindx;
                 Dmat->indx          = d2_indx;    Dmat->val      = d2_inv;
                 Dmat->data_org      = data_org;
                 Dmat->matvec        = precond->Pmat->matvec;
                 Dmat->matrix_type   = precond->Pmat->matrix_type;

                 if (options[AZ_pre_calc] != AZ_reuse) {
                    AZ_calc_blk_diag_inv(val, indx, bindx, rpntr, cpntr, bpntr,
                                         d2_inv, d2_indx, d2_bindx, d2_rpntr, 
                                         d2_bpntr, data_org);
                 }
                 else if (i == AZ_NEW_ADDRESS) {
                   AZ_printf_err( "Error: options[AZ_pre_calc]==AZ_reuse and"
                         "previous factors\n       not found. Check"
                         "data_org[AZ_name].\n");
                   exit(-1);
                 }
           }
           else if (previous_factors != data_org[AZ_name]) {
              AZ_printf_err( "Warning: Using a previous factorization as a"
                       "preconditioner\neven though matrix"
                       "(data_org[AZ_name]) has changed\n");
           }
           previous_factors = data_org[AZ_name];

           /* scale rhs */

           sprintf(tag,"v_prec %s",precond->context->tag);
           v = AZ_manage_memory((N+max_externals)*sizeof(double),
                           AZ_ALLOC, AZ_SYS+az_iterate_id, tag, &i);

           Dmat->matvec(current_rhs, v, Dmat, proc_config);

           DCOPY_F77(&N, v, &ione, current_rhs, &ione);

           if (options[AZ_poly_ord] > 1) {
              sprintf(tag,"y_prec %s",precond->context->tag);
              y = AZ_manage_memory((N+max_externals)*sizeof(double),
                             AZ_ALLOC, AZ_SYS+az_iterate_id, tag, &i);

              sprintf(tag,"temp_prec %s",precond->context->tag);
              temp = AZ_manage_memory(N*sizeof(double), AZ_ALLOC,AZ_SYS+az_iterate_id,tag,&i);

              for (step = 1; step < options[AZ_poly_ord]; step++) {
                 Amat->matvec(v, y, Amat, proc_config);
                 Dmat->matvec(y, temp, Dmat, proc_config);

                 for (i = 0; i < N; i++) v[i] += current_rhs[i] - temp[i];
              }

              for (i = 0; i < N; i++) current_rhs[i] = v[i];
           }
        }
     break;
     case AZ_sym_GS:

        /* symmetric Gauss-Seidel preconditioner only available on 1 proc */

        if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) AZ_sym_gauss_seidel();
        else if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX)
           AZ_sym_gauss_seidel_sl(val, bindx, current_rhs, data_org, options,
				  precond->context, proc_config);
     break;

     case AZ_Neumann:
     case AZ_ls:
        if (!options[AZ_poly_ord]) return;
        AZ_polynomial_expansion(current_rhs, options, proc_config, precond);
     break;

     case AZ_dom_decomp:
     case AZ_rilu:
        AZ_domain_decomp(current_rhs, precond->Pmat, options, proc_config, 
                         params, precond->context);
     break;

     case AZ_icc:
        /* incomplete Cholesky factorization */

        (void) AZ_printf_out("Incomplete Cholesky not available (use ilu).\n");
     break;

     case AZ_user_precond:
        precond->prec_function(current_rhs, options, proc_config, 
                               params, Amat, precond);
     break;
     case AZ_smoother:
        sprintf(label,"istatus %s",precond->context->tag);
        istatus = AZ_manage_memory(AZ_STATUS_SIZE*sizeof(double),AZ_ALLOC,
				   AZ_SYS+az_iterate_id, label,&i);
        for (i = 0 ; i < AZ_STATUS_SIZE ; i++ ) istatus[i] = 0.0;

        sprintf(label,"y %s",precond->context->tag);
        y = AZ_manage_memory((N+max_externals)*sizeof(double), AZ_ALLOC, 
			     AZ_SYS+az_iterate_id, label, &i);
        sprintf(label,"tttemp %s",precond->context->tag);
        tttemp = AZ_manage_memory((N+max_externals)*sizeof(double),AZ_ALLOC,
				  AZ_SYS+az_iterate_id, label, &i);

        for (i = 0 ; i < N ; i++ ) tttemp[i] = current_rhs[i];

        N_fixed = 0; fixed_pts = NULL;
        if (Amat->aux_ival != NULL) {
           N_fixed   = Amat->aux_ival[0][0];
           fixed_pts = Amat->aux_ival[1];
        }
        else if (options[AZ_pre_calc] != AZ_sys_reuse)
           AZ_printf_out("Warning: Not fixed points set for local smoothing!!\n");

        for (j = 0; j < options[AZ_poly_ord]; j++) {
           AZ_loc_avg(Amat, tttemp, y, N_fixed, fixed_pts, proc_config);
           norm1 = sqrt(AZ_gdot(N, y, y, proc_config));
           if (proc_config[AZ_node] == 0) {
              if ((j==0) && (options[AZ_output] != AZ_none) &&
                  (options[AZ_output] != AZ_last) &&
                  (options[AZ_output] != AZ_summary) &&
                  (options[AZ_output] != AZ_warnings))
                  AZ_printf_out("   %d  %e\n",j, norm1);
              else if ((j==options[AZ_poly_ord]-1) && 
		  (options[AZ_output] != AZ_none) && 
                  (options[AZ_output] != AZ_warnings))
                  AZ_printf_out("   %d  %e\n",j, norm1);
              else if ((options[AZ_output] > 0) && (j%options[AZ_output] == 0))
                  AZ_printf_out("   %d  %e\n",j, norm1);
           }
           for (i = 0 ; i < N ; i++ ) tttemp[i] = y[i];
        }
        for (i = 0 ; i < N ; i++ ) y[i] = current_rhs[i] - y[i];
        for (i = 0 ; i < N ; i++ ) current_rhs[i] = 0.0;

        opt_save1 = options[AZ_output];
        opt_save2 = options[AZ_solver];
        opt_save3 = options[AZ_precond];
        opt_save4 = options[AZ_max_iter];
        opt_save5 = options[AZ_aux_vec];

        options[AZ_output]  = AZ_warnings;
        options[AZ_solver]  = AZ_tfqmr;
        options[AZ_precond] = AZ_dom_decomp;
        options[AZ_max_iter]= 1000;
        options[AZ_aux_vec] = AZ_rand;

        options[AZ_recursion_level]++;
        AZ_oldsolve(current_rhs, y,options, params, istatus, proc_config, 
                    Amat, precond, NULL);
        options[AZ_recursion_level]--;
        options[AZ_output]  = opt_save1;
        options[AZ_solver]  = opt_save2;
        options[AZ_precond] = opt_save3;
        options[AZ_max_iter]= opt_save4;
        options[AZ_aux_vec] = opt_save5;
     break;
     default:
        if (options[AZ_precond] < AZ_SOLVER_PARAMS) {
           AZ_recover_sol_params(options[AZ_precond], &ioptions, &iparams,
                                 &istatus, &Aptr, &Pptr, &Sptr);
           sprintf(label,"y %s",precond->context->tag);
           y = AZ_manage_memory((N+max_externals)*sizeof(double),
                                AZ_ALLOC, AZ_SYS+az_iterate_id, label, &i);
           for (i = 0 ; i < N ; i++ ) y[i] = current_rhs[i];
           for (i = 0 ; i < N ; i++ ) current_rhs[i] = 0.0;

           ioptions[AZ_recursion_level] = options[AZ_recursion_level] + 1;
           if ((options[AZ_pre_calc] == AZ_sys_reuse) &&
               (ioptions[AZ_keep_info] == 1)) 
              ioptions[AZ_pre_calc] = AZ_reuse;
           AZ_oldsolve(current_rhs, y,ioptions,iparams, istatus, proc_config, 
                       Aptr, Pptr, Sptr);
        }
        else {
           (void) AZ_printf_err( "%sERROR: invalid preconditioning flag.\n"
                   "       options[AZ_precond] improperly set (%d).\n", yo,
			   options[AZ_precond]);
           exit(-1);
        }

     }
     options[AZ_pre_calc] = AZ_sys_reuse;
     precond->context->Pmat_computed = 1;

     if (multilevel_flag) {
        if (precond->next_prec == NULL) {
           multilevel_flag = 0;
           for (i = 0; i < N; i++) current_rhs[i] += x_precond[i];
        }
        else {
           for (i = 0; i < N; i++) x_precond[i] += current_rhs[i];
           AZ_compute_residual(orig_rhs, x_precond, current_rhs, 
                               proc_config, Amat);
           precond = precond->next_prec;
           options = precond->options;
           params  = precond->params;
        }
     }

  } while (multilevel_flag);

  proc_config[AZ_MPI_Tag] = AZ_MSG_TYPE;   /* reset all the message types.   */
                                           /* This is to make sure that all  */
                                           /* processors (even those without */
                                           /* any preconditioning work) have */
                                           /* the same message types for the */
                                           /* next message.                  */
#ifdef TIMING
  ttt = AZ_second() - ttt;
  if (input_options[AZ_recursion_level] == 0) input_precond->timing[0] += ttt;
#endif

} /* precond */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_calc_blk_diag_inv(double *val, int *indx, int *bindx, int *rpntr,
                          int *cpntr, int *bpntr, double *d_inv, int *d_indx,
                          int *d_bindx, int *d_rpntr, int *d_bpntr,
                          int *data_org)

/*******************************************************************************

  Routine to calculate the inverse of the block-diagonal portion of the sparse
  matrix in 'val' and the associated integer pointer vectors. This is used for
  scaling and/or preconditioning.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                   Aztec User's Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  d_inv:           Vector containing the inverses of the diagonal blocks.

  d_indx:          The 'indx' array corresponding to the inverse-block
                   diagonals.

  d_bindx:         The 'bindx' array corresponding to the inverse-block
                   diagonals.

  d_rpntr:         The 'rpntr' array corresponding to the inverse-block
                   diagonals.

  d_bpntr:         The 'bpntr' array corresponding to the inverse-block
                   diagonals.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

*******************************************************************************/

{

  /* local variables */

  register int i, j, iblk_row, jblk, icount = 0, iblk_count = 0, ival;
  int          m1, n1, itemp;
  int          m;
  int          bpoff, idoff;
  int         *ipiv, info;
  double      *work;
  char        *yo = "AZ_calc_blk_diag_inv: ";

  /**************************** execution begins ******************************/

  m = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

  if (m == 0) return;

  /* allocate vectors for lapack routines */

  ipiv = (int *)    AZ_allocate(rpntr[m]*sizeof(int));
  work = (double *) AZ_allocate(rpntr[m]*sizeof(double));
  if (work == NULL) AZ_perror("Not enough space for Block Jacobi\n");

  /* offset of the first block */

  bpoff = *bpntr;
  idoff = *indx;

  /* loop over block rows */

  for (iblk_row = 0; iblk_row < m; iblk_row++) {

    /* number of rows in the current row block */

    m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

    /* starting index of current row block */

    ival = indx[bpntr[iblk_row] - bpoff] - idoff;

    /* loop over column block numbers, looking for the diagonal block */

    for (j = bpntr[iblk_row] - bpoff; j < bpntr[iblk_row+1] - bpoff; j++) {
      jblk = bindx[j];

      /* determine the number of columns in this block */

      n1 = cpntr[jblk+1] - cpntr[jblk];

      itemp = m1*n1;

      if (jblk == iblk_row) {   /* diagonal block */

        /* error check */

        if (n1 != m1) {
          (void) AZ_printf_err( "%sERROR: diagonal blocks are not square\n.",
                         yo);
          exit(-1);
        }
        else {

          /* fill the vectors */

          d_indx[iblk_count]  = icount;
          d_rpntr[iblk_count] = rpntr[iblk_row];
          d_bpntr[iblk_count] = iblk_row;
          d_bindx[iblk_count] = iblk_row;

          for (i = 0; i < itemp; i++) d_inv[icount++] = val[ival + i];

          /* invert the dense matrix */

          DGETRF_F77(&m1, &m1, &d_inv[d_indx[iblk_count]], &m1, ipiv, &info);

          if (info < 0) {
            (void) AZ_printf_err( "%sERROR: argument %d is illegal.\n", yo,
                           -info);
            exit(-1);
          }

          else if (info > 0) {
            (void) AZ_printf_err( "%sERROR: the factorization has produced a "
                           "singular U with U[%d][%d] being exactly zero.\n",
                           yo, info, info);
            exit(-1);
          }

          DGETRI_F77(&m1, &d_inv[d_indx[iblk_count]], &m1, ipiv, work, &m1, &info);

          if (info < 0) {
            (void) AZ_printf_err( "%sERROR: argument %d is illegal.\n", yo,
                           -info);
            exit(-1);
          }

          else if (info > 0) {
            (void) AZ_printf_err( "%sERROR: U[%d][%d] is exactly zero;\n", yo,
                           info, info);
            (void) AZ_printf_err( "the matrix is singular and its inverse "
                           "could not be computed.\n");
            exit(-1);
          }
          iblk_count++;
        }
        break;
      }
      else
        ival += itemp;
    }
  }

  d_indx[iblk_count]  = icount;
  d_rpntr[iblk_count] = rpntr[iblk_row];
  d_bpntr[iblk_count] = iblk_row;

  AZ_free((void *) ipiv);
  AZ_free((void *) work);

} /* AZ_calc_blk_diag_inv */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_calc_blk_diag_LU(double *val, int *indx, int *bindx, int *rpntr,
                          int *cpntr, int *bpntr, double *d_inv, int *d_indx,
                          int *d_bindx, int *d_rpntr, int *d_bpntr,
                          int *data_org, int *ipvt)

/*******************************************************************************

  Routine to calculate the LU factors of the block-diagonal portion of sparse
  matrix in 'val' and the associated integer pointer vectors. This is used for
  scaling.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                   Aztec User's Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  d_inv:           Vector containing the LU of the diagonal blocks.

  d_indx:          The 'indx' array corresponding to the LU-block
                   diagonals.

  d_bindx:         The 'bindx' array corresponding to the LU-block
                   diagonals.

  d_rpntr:         The 'rpntr' array corresponding to the LU-block
                   diagonals.

  d_bpntr:         The 'bpntr' array corresponding to the LU-block
                   diagonals.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

*******************************************************************************/

{

  /* local variables */

  register int i, j, iblk_row, jblk, icount = 0, iblk_count = 0, ival;
  int          m1, n1, itemp;
  int          m;
  int          bpoff, idoff;
  int         info;
  double      *work;
  char        *yo = "AZ_calc_blk_diag_inv: ";

  /**************************** execution begins ******************************/

  m = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

  if (m == 0) return;

  /* allocate vectors for lapack routines */

  work = (double *) AZ_allocate(rpntr[m]*sizeof(double));
  if (work == NULL) AZ_perror("Not enough space for Block Jacobi\n");

  /* offset of the first block */

  bpoff = *bpntr;
  idoff = *indx;

  /* loop over block rows */

  for (iblk_row = 0; iblk_row < m; iblk_row++) {

    /* number of rows in the current row block */

    m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

    /* starting index of current row block */

    ival = indx[bpntr[iblk_row] - bpoff] - idoff;

    /* loop over column block numbers, looking for the diagonal block */

    for (j = bpntr[iblk_row] - bpoff; j < bpntr[iblk_row+1] - bpoff; j++) {
      jblk = bindx[j];

      /* determine the number of columns in this block */

      n1 = cpntr[jblk+1] - cpntr[jblk];

      itemp = m1*n1;

      if (jblk == iblk_row) {   /* diagonal block */

        /* error check */

        if (n1 != m1) {
          (void) AZ_printf_err( "%sERROR: diagonal blocks are not square\n.",
                         yo);
          exit(-1);
        }
        else {

          /* fill the vectors */

          d_indx[iblk_count]  = icount;
          d_rpntr[iblk_count] = rpntr[iblk_row];
          d_bpntr[iblk_count] = iblk_row;
          d_bindx[iblk_count] = iblk_row;

          for (i = 0; i < itemp; i++) d_inv[icount++] = val[ival + i];

          /* invert the dense matrix */

          DGETRF_F77(&m1, &m1, &d_inv[d_indx[iblk_count]], &m1, &(ipvt[rpntr[iblk_row]]), &info);

          if (info < 0) {
            (void) AZ_printf_err( "%sERROR: argument %d is illegal.\n", yo,
                           -info);
            exit(-1);
          }

          else if (info > 0) {
            (void) AZ_printf_err( "%sERROR: the factorization has produced a "
                           "singular U with U[%d][%d] being exactly zero.\n",
                           yo, info, info);
            exit(-1);
          }
          iblk_count++;
        }
        break;
      }
      else
        ival += itemp;
    }
  }

  d_indx[iblk_count]  = icount;
  d_rpntr[iblk_count] = rpntr[iblk_row];
  d_bpntr[iblk_count] = iblk_row;

  AZ_free((void *) work);

} /* AZ_calc_blk_diag_inv */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void jacobi(double val[], double x[], int data_org[])

/*******************************************************************************

  Simple Jacobi iteration (undamped).  Not yet implemented for DVBR formatted
  matrices.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                   Aztec User's Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  x:               On input, contains the current solution to the linear system.
                   On output contains the Jacobi preconditioned solution.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

*******************************************************************************/

{

  /* local variables */

  register int i;
  int          N;

  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  for (i = 0; i < N; i++)
    *x++ /=  *val++;

} /* jacobi */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

extern void AZ_sym_gauss_seidel(void)

/*******************************************************************************

  Symmetric Gauss-Siedel preconditioner

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:  void
  ===============

*******************************************************************************/

{

  /* local variables */

  /**************************** execution begins ******************************/

  (void) AZ_printf_err( "WARNING: sym Gauss-Seidel preconditioning not\n"
                 "         implemented for VBR matrices\n");
  exit(-1);

} /* AZ_sym_gauss_seidel */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sym_gauss_seidel_sl(double val[],int bindx[],double x[],int data_org[],
			    int options[], struct context *context,
			    int proc_config[])

/*******************************************************************************

  Symmetric Gauss-Siedel preconditioner.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                   Aztec User's Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  x:               On input, contains the current solution to the linear system.
                   On output contains the Jacobi preconditioned solution.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

  options:         Determines specific solution method and other parameters.

*******************************************************************************/

{

  /* local variables */

  register int    *bindx_ptr;
  register double sum, *ptr_val;
  int             i, bindx_row, j_last, N, step, ione = 1, j;
  double          *b, *ptr_b;
  char            tag[80];

  /**************************** execution begins ******************************/

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  sprintf(tag,"b/sGS %s",context->tag);
  b = AZ_manage_memory(N*sizeof(double), AZ_ALLOC, AZ_SYS+az_iterate_id, tag, &i);

  DCOPY_F77(&N, x, &ione, b, &ione);
  ptr_val = val;

  for (i = 0; i < N; i++) {
    (*ptr_val) = 1.0 / (*ptr_val);
    x[i]     = 0.0;
    ptr_val++;
  }

  for (step = 0; step < options[AZ_poly_ord]; step++) {
    AZ_exchange_bdry(x, data_org, proc_config);

    bindx_row = bindx[0];
    bindx_ptr = &bindx[bindx_row];
    ptr_val   = &val[bindx_row];
    ptr_b   = b;

    for (i = 0; i < N; i++) {
      sum    = *ptr_b++;
      j_last = bindx[i+1] - bindx[i];

      for (j = 0; j < j_last; j++) {
        sum -= *ptr_val++ * x[*bindx_ptr++];
      }
      x[i] = sum * val[i];
    }

    bindx_row = bindx[N];
    bindx_ptr = &bindx[bindx_row-1];
    ptr_val   = &val[bindx_row-1];

    for (i = N - 1; i >= 0; i--) {
      sum = b[i];
      j_last  = bindx[i+1] - bindx[i];

      for (j = 0; j < j_last; j++) {
        sum -= *ptr_val-- * x[*bindx_ptr--];
      }
      x[i] = sum * val[i];
    }
  }

  for (i = 0; i < N; i++)
    val[i] = 1.0 / val[i];

} /* AZ_sym_gauss_seidel_sl */

#ifdef eigen
void AZ_do_Jacobi(double val[], int indx[], int bindx[], int rpntr[],
                     int cpntr[], int bpntr[], double x[], double b[],
                     double temp[], int options[], int data_org[],
                     int proc_config[], double params[], int flag)
{
double *v;
int i,step;
int N;
 
  N    = data_org[AZ_N_internal] + data_org[AZ_N_border];
 
    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
 
      if ( (options[AZ_poly_ord] != 0) && (flag == 1) )
         for (i = data_org[AZ_N_internal]; i < N; i++) x[i] = b[i]/val[i];
 
      if (options[AZ_poly_ord] > flag) {
        v = AZ_manage_memory((N+data_org[AZ_N_external])*sizeof(double),
                             AZ_ALLOC, AZ_SYS+az_iterate_id, "v in do_jacobi", &i);
 
        for (i = 0; i < N; i++) v[i] = x[i];
 
        for (step = flag; step < options[AZ_poly_ord]; step++) {
          Amat->matvec(v, temp, Amat, proc_config);

          for(i = 0; i < N; i++) v[i] += (b[i] - temp[i]) / val[i];
        }
        for (i = 0; i < N; i++) x[i] = v[i];
      }
    }
    else {
       (void) AZ_printf_err("AZ_slu with option[AZ_poly_ord] > 0 only \n");
       (void) AZ_printf_err("implemented for MSR matrices.\n");
       exit(-1);
    }

}
#endif
