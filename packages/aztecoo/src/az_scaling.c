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
#include <string.h>
#include "az_aztec.h"
#include "az_blas_wrappers.h"
#include "az_lapack_wrappers.h"

/* static functions */

static void calc_blk_diag_Chol(double *val, int *indx, int *bindx, int *rpntr,
                               int *cpntr, int *bpntr, double *L, int *d_indx,
                               int *d_bindx, int *d_rpntr, int *d_bpntr,
                               int *data_org);

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_scale_f(int action, AZ_MATRIX *Amat, int options[], double b[], double x[], 
	        int proc_config[], struct AZ_SCALING *scaling)

/*******************************************************************************

  Scale the matrix, rhs and initial guess as specified by the user.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                    Aztec Users Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file  Aztec Users Guide).

  b:               Right hand side of linear system.

  x:               On input, contains the initial guess. On output contains the
                   solution to the linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see  Aztec Users Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  action:          Flag which determines whether to scale the matrix/rhs/sol or 
                   just the rhs or just the solution or to inverse scale the
                   rhs or solution.

  Amat:            Structure used to represent the matrix (see az_aztec.h
                   and Aztec User's Guide).


*******************************************************************************/

{

  /* local variables */

  char *yo = "AZ_scale: ";
  int  temp;

  /**************************** execution begins ******************************/

  temp = options[AZ_scaling];
  if ( (options[AZ_scaling] != AZ_none) && 
       (Amat->data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
       (Amat->data_org[AZ_matrix_type] != AZ_VBR_MATRIX)  ) {
    if (scaling->scale != 0) {
      int err = scaling->scale(action, Amat, options, b, x, proc_config, scaling);
      if (err != 0) {
        AZ_printf_out("AZ_scale_f ERROR, code %d returned from user-scaling function\n",
               err);
      }
      return;
    }
    else {
      AZ_printf_out("AZ_scale_f WARNING: have user-matrix, but scaling struct contains null 'scale' callback. Turning off scaling.\n");
      options[AZ_scaling] = AZ_none;
    }
  }

  if ( (options[AZ_ignore_scaling] == AZ_TRUE) &&
       (options[AZ_pre_calc] <= AZ_recalc) &&
       (action == AZ_SCALE_MAT_RHS_SOL)) {
    scaling->A_norm = AZ_gmax_matrix_norm(Amat->val, Amat->indx, Amat->bindx, 
                                         Amat->rpntr, Amat->cpntr, Amat->bpntr, 
                                          proc_config, Amat->data_org);
  }

  switch (options[AZ_scaling]) {

  case AZ_none:
    break;

  case AZ_Jacobi:
    AZ_block_diagonal_scaling(action, Amat, Amat->val, Amat->indx, Amat->bindx, 
                              Amat->rpntr, Amat->cpntr, Amat->bpntr, Amat->data_org, 
			      b, options, proc_config, scaling);
    break;

  case AZ_BJacobi:
    AZ_block_diagonal_scaling(action, Amat, Amat->val, Amat->indx, Amat->bindx, 
                              Amat->rpntr, Amat->cpntr, Amat->bpntr, Amat->data_org, 
			      b, options, proc_config, scaling);
    break;

  case AZ_row_sum:
    AZ_row_sum_scaling(action, Amat, b, options, scaling);
    break;

  case AZ_sym_diag:
    AZ_sym_diagonal_scaling(action,Amat,b,x,options,proc_config,scaling);
    break;

  case AZ_sym_row_sum:
    AZ_sym_row_sum_scaling(action, Amat, 
                              b, x, options, proc_config, scaling);
    break;

  case AZ_equil:
    AZ_equil_scaling(action, Amat,
                              b, x, options, proc_config, scaling);
    break;

  case AZ_sym_BJacobi:

    /* symmetric block Jacobi scaling */
#ifdef next_version
    AZ_sym_block_diagonal_scaling(val, indx, bindx, rpntr, cpntr, bpntr, b,
                                    options, data_org, proc_config);
    else if (action == AZ_RESCALE_SOL)
      AZ_sym_rescale_vbr(x, data_org);
#endif
    break;

  default:
    (void) AZ_printf_err( "%sERROR: invalid scaling option.\n"
                   "          options[AZ_scaling] = %d\n", yo,
                   options[AZ_scaling]);
  exit(-1);
  }
  options[AZ_scaling] = temp;

} /* AZ_scale */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_block_diagonal_scaling(int action, AZ_MATRIX *Amat, double val[], int indx[], 
			       int bindx[], int rpntr[], int cpntr[], int bpntr[], 
			       int data_org[], double b[], int options[], int proc_config[],
			       struct AZ_SCALING *scaling)

/*******************************************************************************

  Routine to block Jacobi scale the sparse matrix problem.  Note: this scales
  the matrix and the right-hand side and the resulting matrix is non-symmetric.

  If the matrix is in MSR format, it is treated as point entry (block size 1)
  and standard Jacobi scaling is performed.  Else, if the matrix is in VBR
  format, block Jacobi scaling is performed.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                    Aztec Users Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file  Aztec Users Guide).

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see  Aztec Users Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  action:          Flag which determines whether to scale the matrix or to
                   rescale the solution.

  Amat:            Structure used to represent the matrix (see az_aztec.h
                   and Aztec User's Guide).


*******************************************************************************/

{

  /* local variables */

  static AZ_MATRIX *Amat_inv;
  register int   iblk_row, i, j, k, ival, d_ival, jblk, ib;
  int            bpoff, idoff;
  int            d_bpoff, d_idoff;
  int            m1, n1, ib1, ib2;
  int            ione = 1, itemp;
  int            m, Proc, N;
  int            tsize;
  int            j_last, bindx_row;
  int            max_blk;
  static int    *d3_indx, *d3_bindx, *d3_rpntr, *d3_bpntr;
  static double *d3_inv, *sc_vec;
  double         *c, *work;
  char           None[2], Left[2], Lower[2], Unit[2], Upper[2];
  char          *yo = "AZ_block_diagonal_scaling: ";
  char          label[80];
int *ipiv,info, iminus_one = -1;
double one = 1.0;

  /**************************** execution begins ******************************/

  if ( (action == AZ_INVSCALE_SOL) || (action == AZ_SCALE_SOL)) return;

  /* initialize */

  N    = data_org[AZ_N_internal] + data_org[AZ_N_border];
  m    = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];
  Proc = proc_config[AZ_node];

  scaling->action = AZ_left_scaling;

  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {

    /***** MSR Matrix => point Jacobi scaling *****/

    sprintf(label,"sc_vec%d",options[AZ_recursion_level]);
    sc_vec = (double *) AZ_manage_memory((N+data_org[AZ_N_external])*sizeof(double), AZ_ALLOC,
                                         data_org[AZ_name], label,
                                         &itemp);

    if (((action == AZ_SCALE_RHS) || (action == AZ_INVSCALE_RHS) || 
         (options[AZ_pre_calc] >= AZ_reuse)) && (itemp == AZ_NEW_ADDRESS)) {

     if (itemp == AZ_NEW_ADDRESS) {
        AZ_manage_memory(0, AZ_CLEAR, data_org[AZ_name], label, (int*)0);

        /* If we're about to exit due to getting itemp==AZ_NEW_ADDRESS, let's try
           one other thing. See if the required memory can be found using
           type==scaling->mat_name. If not, then go ahead and crash. ABW
        */
        if (scaling != 0) {
          sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                               sizeof(double), AZ_ALLOC,
                                               scaling->mat_name, label, &itemp);
        }
      }

      if (itemp == AZ_NEW_ADDRESS) {
        (void) AZ_printf_err( "%sERROR: Previous scaling not found for matrix "
                     "with\ndata_org[AZ_name] = %d. Check\n"
                     "options[AZ_pre_calc]\n", yo,
                     data_org[AZ_name]);
        exit(-1);
      }
    }

    if ( (options[AZ_pre_calc] <= AZ_recalc) && (action == AZ_SCALE_MAT_RHS_SOL)) {
      for (ib = 0; ib < N; ib++) {
        if (fabs(val[ib]) > DBL_MIN) sc_vec[ib] = 1.0 / val[ib];
        else                         sc_vec[ib] = 1.0;

        val[ib] = 1.0;
        j_last  = bindx[ib+1] - bindx[ib];
        bindx_row = bindx[ib];

        for (j = 0; j < j_last; j++) {
          k       = bindx_row + j;
          val[k] *= sc_vec[ib];
        }
      }
    }
    if ( (action == AZ_SCALE_MAT_RHS_SOL) || (action == AZ_SCALE_RHS) ) {
       for (ib = 0; ib < N; ib++) b[ib] *= sc_vec[ib];
    }
    if ( action == AZ_INVSCALE_RHS)  {
       for (ib = 0; ib < N; ib++) b[ib] /= sc_vec[ib];
    }
  }

  else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {

    /***** VBR Matrix => block Jacobi scaling *****/

    sprintf(None, "N");

    /* First, compute the block-diagonal inverse (if not already computed) */

    tsize = 0;
    for (i = 0; i < m; i++)
      tsize += (rpntr[i+1] - rpntr[i]) * (cpntr[i+1] - cpntr[i]);

    sprintf(label,"d3_indx%d",options[AZ_recursion_level]);
    d3_indx  = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                           data_org[AZ_name], label, &itemp);
    sprintf(label,"d3_bindx%d",options[AZ_recursion_level]);
    d3_bindx = (int *)    AZ_manage_memory(m*sizeof(int), AZ_ALLOC,
                                           data_org[AZ_name], label, &itemp);
    sprintf(label,"d3_rpntr%d",options[AZ_recursion_level]);
    d3_rpntr = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                           data_org[AZ_name], label, &itemp);
    sprintf(label,"d3_bpntr%d",options[AZ_recursion_level]);
    d3_bpntr = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                           data_org[AZ_name], label, &itemp);
    sprintf(label,"d3_inv%d",options[AZ_recursion_level]);
    d3_inv   = (double *) AZ_manage_memory(tsize*sizeof(double), AZ_ALLOC,
                                           data_org[AZ_name], label, &itemp);
    sprintf(label,"Amat_inv%d",options[AZ_recursion_level]);
    Amat_inv = (AZ_MATRIX *) AZ_manage_memory(sizeof(AZ_MATRIX), AZ_ALLOC,
                                           data_org[AZ_name], label, &itemp);

    sprintf(label,"ipv %d",options[AZ_recursion_level]);
    ipiv = (int *) AZ_manage_memory((N+1)*sizeof(int), AZ_ALLOC,
                                           data_org[AZ_name], label, &itemp);

    Amat_inv->rpntr       = d3_rpntr;   Amat_inv->cpntr    = d3_rpntr;
    Amat_inv->bpntr       = d3_bpntr;   Amat_inv->bindx    = d3_bindx;
    Amat_inv->indx        = d3_indx;    Amat_inv->val      = d3_inv;
    Amat_inv->data_org    = data_org;
    Amat_inv->matvec      = Amat->matvec;
    Amat_inv->matrix_type = Amat->matrix_type;


    if (((action == AZ_SCALE_RHS) || (action == AZ_INVSCALE_RHS) || 
         (options[AZ_pre_calc] >= AZ_reuse)) && (itemp == AZ_NEW_ADDRESS)) {
      (void) AZ_printf_err( "%sERROR: Previous scaling not found for matrix "
                     "with\ndata_org[AZ_name] = %d. Check\n"
                     "options[AZ_pre_calc]\n", yo,
                     data_org[AZ_name]);
      exit(-1);
    }



    if ( (options[AZ_pre_calc] <= AZ_recalc) && (action == AZ_SCALE_MAT_RHS_SOL)) {

      AZ_calc_blk_diag_LU(val, indx, bindx, rpntr, cpntr, bpntr, d3_inv,
                           d3_indx, d3_bindx, d3_rpntr, d3_bpntr, data_org, ipiv);

      /* offset of the first block */

      bpoff = *bpntr;
      idoff = *indx;
      d_bpoff = *d3_bpntr;
      d_idoff = *d3_indx;

      /* scale the matrix 'A' */

      max_blk = 0;
      for (i = 0; i < m + data_org[AZ_N_ext_blk] ; i++) {
        if ( cpntr[i+1]-cpntr[i] > max_blk )
          max_blk = cpntr[i+1] - cpntr[i];
      }

      work = (double *) AZ_allocate(max_blk*max_blk*sizeof(double));
      if (work == NULL) {
        (void) AZ_printf_err( "%sERROR: not enough memory for diagonal\n"
                       "      scaling. Not able to allocate work\n"
                       "      array. Must run a smaller problem\n", yo);
        exit(-1);
      }

      /* loop over the block rows */

      for (iblk_row = 0; iblk_row < m; iblk_row++) {

        /* find out how many rows are in this block row */

        m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

        /* starting index of current row block */

        ival = indx[bpntr[iblk_row] - bpoff] - idoff;

        /* starting index of current block row for diagonal scaling blocks */

        d_ival  = d3_indx[d3_bpntr[iblk_row] - d_bpoff] - d_idoff;

        /* loop over the (block) columns in the current (block) row */

        for (j = bpntr[iblk_row] - bpoff; j < bpntr[iblk_row+1] - bpoff; j++){
          jblk = bindx[j];

          /* the starting point column index of the current block */

          ib1 = cpntr[jblk];

          /* ending point column index of the current block */

          ib2 = cpntr[jblk+1];

          /* number of columns in the current block */

          n1    = ib2 - ib1;
          itemp = m1*n1;

          if (jblk == iblk_row) {

            /* diagonal block => set to identity */

            if (m1 != n1) {
              if (Proc == 0) {
                (void) AZ_printf_err( "%sERROR: diagonal blocks are not\n"
                               "       square inside scaling\n", yo);
              }
              exit(-1);
            }

            for (i = 0; i < m1; i++)
              for (k = 0; k < n1; k++)
                if (i == k)
                  val[ival + i + m1*k] = 1.0;
                else
                  val[ival + i + m1*k] = 0.0;
          }

          else {
            if (itemp > max_blk*max_blk) {
              (void) AZ_printf_err( "%sERROR: block size (%d) is too big =>\n",
                             yo, itemp);
              exit(-1);
            }

            /* Matrix solve */

            DGETRS_F77(CHAR_MACRO(None[0]), &m1, &n1, &d3_inv[d_ival], &m1, 
                    &(ipiv[rpntr[iblk_row]]), &val[ival], &m1, &info);
          }
          ival += itemp;
        }
      }

      AZ_free((void *) work);



    }
    d_bpoff = *d3_bpntr;
    d_idoff = *d3_indx;

    /* lastly, scale the rhs */

    c = (double *) AZ_allocate((unsigned) N * sizeof(double));
    if (c == NULL) {
       (void) AZ_printf_err( "%sERROR: not enough memory for block diagonal\n"
                     "       scaling. Not able to allocate c\n"
                     "       array. Must run a smaller problem\n", yo);
       exit(-1);
    }
    if ( action == AZ_INVSCALE_RHS)  {
       sprintf(Left,  "L");
       sprintf(Upper, "U");
       sprintf(Lower, "L");
       sprintf(Unit, "U");
       for (iblk_row = 0; iblk_row < m; iblk_row++) {
          m1 = rpntr[iblk_row+1] - rpntr[iblk_row];
          d_ival  = d3_indx[d3_bpntr[iblk_row] - d_bpoff] - d_idoff;
          DTRMM_F77( CHAR_MACRO(Left[0]), CHAR_MACRO(Upper[0]), CHAR_MACRO(None[0]), CHAR_MACRO(None[0]), &m1,
                  &ione, &one, &d3_inv[d_ival], &m1, &(b[rpntr[iblk_row]]), &m1);
          DTRMM_F77( CHAR_MACRO(Left[0]), CHAR_MACRO(Lower[0]), CHAR_MACRO(None[0]), CHAR_MACRO(Unit[0]), &m1, &ione,
                  &one, &d3_inv[d_ival], &m1, &(b[rpntr[iblk_row]]), &m1);
          AZ_DLASWP_F77( &ione, &(b[rpntr[iblk_row]]), &m1, &ione, &m1, &(ipiv[rpntr[iblk_row]]), &iminus_one );

       }
    }
    if ( (action == AZ_SCALE_MAT_RHS_SOL) || (action == AZ_SCALE_RHS) ) {
       for (iblk_row = 0; iblk_row < m; iblk_row++) {
          m1 = rpntr[iblk_row+1] - rpntr[iblk_row];
          d_ival  = d3_indx[d3_bpntr[iblk_row] - d_bpoff] - d_idoff;
          DGETRS_F77(CHAR_MACRO(None[0]), &m1, &ione, &d3_inv[d_ival], &m1, &(ipiv[rpntr[iblk_row]]), 
                  &(b[rpntr[iblk_row]]), &m1, &info);
       }
    }
    AZ_free((void *) c);

  }

} /* AZ_block_diagonal_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sym_block_diagonal_scaling(double val[], int indx[], int bindx[],
                                   int rpntr[], int cpntr[], int bpntr[],
                                   double b[], int options[], 
                                   int data_org[], int proc_config[]
                                   /* , struct AZ_SCALING *scaling*/)

/*******************************************************************************

  Routine to symmetric block Jacobi scale the sparse matrix problem.  Note: this
  scales the matrix, the solution and the right-hand side.

  If the matrix is in MSR format, it is treated as point entry (block size 1)
  and standard Jacobi scaling is performed.  Else, if the matrix is in VBR
  format, block symmetric Jacobi scaling is performed.

  For the DVBR format, the following system is created from Ax = b:

       (trans(inv(L)) A inv(L) (L x) = trans(inv(L)) b

  where L is the Cholesky factor of the block diagonal portion of A.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                    Aztec Users Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file  Aztec Users Guide).

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see  Aztec Users Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  action:          Flag which determines whether to scale the matrix or to
                   rescale the solution.

*******************************************************************************/

{

  /* local variables */

  register int   iblk_row, i, j, k, ival, jblk;
  register int   it, jt, icount, iwork, ilow, ib = 0;
  int            bpoff, idoff;
  int            d_bpoff, d_idoff;
  int            m1, n1, ib1, ib2, idL;
  int            ione = 1, itemp, itemp2;
  int            m, Proc;
  int            tsize;
  int            max_blk;
  static int    *d3_indx, *d3_bindx, *d3_rpntr, *d3_bpntr;
  static double *L;
  double         done = 1.0;
  double        *work;
  char           None[2];
  char          *side = "L", *uplo = "L", *transa = "N", *diag = "N";
  char          *yo = "sym_AZ_block_diagonal_scaling: ";
  char           label[80];

  /**************************** execution begins ******************************/

  AZ_printf_out("Error: Symmetric block scaling not implemented\n");
  exit(-1);

  /* initialize */

  m    = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];
  Proc = proc_config[AZ_node];

  sprintf(None, "N");

  /*
   * First, compute the block-diagonal Cholesky factorization, its inverse and
   * its inverse's transpose (if not already computed).
   */

  tsize = 0;
  for (i = 0; i < m; i++)
    tsize += (rpntr[i+1] - rpntr[i]) * (cpntr[i+1] - cpntr[i]);

  sprintf(label,"d3_indx%d",options[AZ_recursion_level]);
  d3_indx  = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                         data_org[AZ_name], label, &itemp);
  sprintf(label,"d3_bindx%d",options[AZ_recursion_level]);
  d3_bindx = (int *)    AZ_manage_memory(m*sizeof(int), AZ_ALLOC,
                                         data_org[AZ_name], label, &itemp);
  sprintf(label,"d3_rpntr%d",options[AZ_recursion_level]);
  d3_rpntr = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                         data_org[AZ_name], label, &itemp);
  sprintf(label,"d3_bpntr%d",options[AZ_recursion_level]);
  d3_bpntr = (int *)    AZ_manage_memory((m+1)*sizeof(int), AZ_ALLOC,
                                         data_org[AZ_name],
                                         label, &itemp);
  sprintf(label,"L%d",options[AZ_recursion_level]);
  L        = (double *) AZ_manage_memory(tsize*sizeof(double), AZ_ALLOC,
                                         data_org[AZ_name], label, &itemp);

  if ((options[AZ_pre_calc] >= AZ_reuse) && (itemp == AZ_NEW_ADDRESS)) {
    (void) AZ_printf_err( "%sERROR: Previous scaling not found for matrix "
                   "with\ndata_org[AZ_name] = %d. Check\n"
                   "options[AZ_pre_calc]\n", yo, data_org[AZ_name]);
    exit(-1);
  }

  /*
   * If necessary, calculate the Cholesky factors (L) of the diagonal blocks and
   * store in L and the d3_ pointer vectors.
   */

  if (options[AZ_pre_calc] <= AZ_recalc) {
    calc_blk_diag_Chol(val, indx, bindx, rpntr, cpntr, bpntr, L, d3_indx,
                       d3_bindx, d3_rpntr, d3_bpntr, data_org);

    /* offset of the first block */

    bpoff = *bpntr;
    idoff = *indx;

    d_bpoff = *d3_bpntr;
    d_idoff = *d3_indx;

    /* symmetrically scale the matrix 'A' using the Cholesky factors */

    max_blk = 0;
    for (i = 0; i < m + data_org[AZ_N_ext_blk] ; i++) {
      if ( cpntr[i+1]-cpntr[i] > max_blk )
        max_blk = cpntr[i+1] - cpntr[i];
    }
    work = (double *) AZ_allocate(max_blk*max_blk*sizeof(double));
    if (work == NULL) {
      (void) AZ_printf_err( "%sERROR: not enough memory for diagonal\n"
                     "      scaling. Not able to allocate work\n"
                     "      array. Must run a smaller problem\n", yo);
      exit(-1);
    }

    /* loop over the block rows */

    for (iblk_row = 0; iblk_row < m; iblk_row++) {

      /* find out how many rows are in this block row */

      m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

      /* starting index of current row block */

      ival = indx[bpntr[iblk_row] - bpoff] - idoff;

      /* loop over the (block) columns in the current (block) row */

      for (j = bpntr[iblk_row] - bpoff; j < bpntr[iblk_row+1] - bpoff; j++){
        jblk = bindx[j];

        /* the starting point column index of the current block */

        ib1 = cpntr[jblk];

        /* ending point column index of the current block */

        ib2 = cpntr[jblk+1];

        /* number of columns in the current block */

        n1    = ib2 - ib1;
        itemp = m1*n1;

        if (jblk == iblk_row) {

          /* diagonal block => set to identity */

          if (m1 != n1) {
            if (Proc == 0) {
              (void) AZ_printf_err( "%sERROR: diagonal blocks are not\n"
                             "       square inside scaling\n", yo);
            }
            exit(-1);
          }

          for (i = 0; i < m1; i++) {
            itemp2 = ival + i;
            for (k = 0; k < n1; k++)
              if (i == k)
                val[itemp2 + m1*k] = 1.0;
              else
                val[itemp2 + m1*k] = 0.0;
          }
        }

        else if (jblk < iblk_row) {

          /* lower off-diagonal block */

          if (itemp > max_blk*max_blk) {
            (void) AZ_printf_err( "%sERROR: block size (%d) is too big =>\n",
                           yo, itemp);
            exit(-1);
          }

          /*
           * Fill the work array with the proper block in "A" and take its
           * transpose.
           */

          icount = 0;
          for (it = 0; it < m1; it++) {
            for (jt = 0; jt < n1; jt++) {
              work[icount] = val[ival+icount];
              icount++;
            }
          }

          AZ_dtrans(&m1, &n1, work);

          /* starting index for the diagonal L block for this first operation */

          idL = d3_indx[d3_bpntr[jblk] - d_bpoff] - d_idoff;

          /* perform a backsolve on L*work' = A' to get 'work' array */

          DTRSM_F77(CHAR_MACRO(side[0]), CHAR_MACRO(uplo[0]), 
		    CHAR_MACRO(transa[0]), CHAR_MACRO(diag[0]), &m1, &n1, 
		    &done, L+idL, &m1, work, &m1);

          /* need the transpose of the work array */

          AZ_dtrans(&m1, &n1, work);

          /* starting index for the diagonal 'L' block for the next operation */

          idL = d3_indx[d3_bpntr[iblk_row] - d_bpoff] - d_idoff;

          /* perform a backsolve on L*work2 = work */

          DTRSM_F77(CHAR_MACRO(side[0]), CHAR_MACRO(uplo[0]), 
		    CHAR_MACRO(transa[0]), CHAR_MACRO(diag[0]), &m1, &n1,
		    &done, L+idL, &m1, work, &m1);

          /* copy the transpose of this result into the proper block in 'val' */

          iwork  = AZ_get_sym_indx(iblk_row, jblk, indx, bindx, bpntr);
          icount = 0;
          for (i = 0; i < n1; i++)
            for (k = 0; k < m1; k++)
              *(val + iwork + i + k*n1) = *(work + icount++);
        }

        ival += itemp;
      }

      /* scale the RHS */

      idL = d3_indx[d3_bpntr[iblk_row] - d_bpoff] - d_idoff;
      DTRSM_F77(CHAR_MACRO(side[0]), CHAR_MACRO(uplo[0]), 
		CHAR_MACRO(transa[0]), CHAR_MACRO(diag[0]), &m1, &ione, &done,
		L+idL, &m1, b+ib, &m1);
      ib += m1;
    }

    AZ_free((void *) work);

    /*
     * Copy the lower block to their corresponding upper blocks for symmetry.
     */

    /* loop over the block rows */

    for (iblk_row = 0; iblk_row < m; iblk_row++) {

      /* find out how many rows are in this block row */

      m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

      /* starting index of current row block */

      ival = indx[bpntr[iblk_row] - bpoff] - idoff;

      /* loop over the (block) columns in the current (block) row */

      for (j = bpntr[iblk_row] - bpoff; j < bpntr[iblk_row+1] - bpoff; j++){
        jblk = bindx[j];

        /* the starting point column index of the current block */

        ib1 = cpntr[jblk];

        /* ending point column index of the current block */

        ib2 = cpntr[jblk+1];

        /* number of columns in the current block */

        n1    = ib2 - ib1;
        itemp = m1*n1;

        if (jblk > iblk_row) {

          /*
           * Lower off-diagonal block => copy to corresponding upper block.
           * NOTE: we do not pass in cpntr as the matrix should be symmetric.
           */

          ilow = AZ_get_sym_indx(iblk_row, jblk, indx, bindx, bpntr);

          icount = 0;
          for (i = 0; i < n1; i++)
            for (k = 0; k < m1; k++)
              val[ilow + i + k*n1] = val[ival + icount++];
        }

        ival += itemp;
      }
    }
  }

} /* sym_AZ_block_diagonal_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_row_sum_scaling(int action, AZ_MATRIX *Amat, double b[], 
                        int options[], struct AZ_SCALING *scaling)

/*******************************************************************************

  Routine to left row-sum scale the sparse matrix problem; Note: this scales
  the entire matrix problem Ax = b

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                    Aztec Users Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file  Aztec Users Guide).

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see  Aztec Users Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  register int indx_ptr, irow, iblk_row, jcol, jblk_col, icol_blk, iblk = 0;
  int          num_blk_rows, num_col_blks, num_blk_cols, I, J, t_index;
  int          m, k, N;
  int          j_last, bindx_row;
  double       row_sum = 0.0, sign = 0.0, inv_row_sum, *sc_vec;
  char        *yo = "AZ_row_sum_scaling: ";
  char        label[80];
  int         *indx, *bindx, *rpntr, *cpntr, *bpntr, *data_org;
  double      *val;

  /**************************** execution begins ******************************/

  if ( (action == AZ_INVSCALE_SOL) || (action == AZ_SCALE_SOL)) return;

  val  = Amat->val;
  indx = Amat->indx;
  bindx = Amat->bindx;
  rpntr = Amat->rpntr;
  cpntr = Amat->cpntr;
  bpntr = Amat->bpntr;
  data_org = Amat->data_org;
  N      = data_org[AZ_N_internal] + data_org[AZ_N_border];

  sprintf(label,"sc_vec%d",options[AZ_recursion_level]);
  sc_vec = (double *) AZ_manage_memory((N+data_org[AZ_N_external])*sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], label, &k);
  scaling->action = AZ_left_scaling;

    if (((action == AZ_SCALE_RHS) || (action == AZ_INVSCALE_RHS) || 
         (options[AZ_pre_calc] >= AZ_reuse)) && (k == AZ_NEW_ADDRESS)) {

      if (k == AZ_NEW_ADDRESS) {
        AZ_manage_memory(0, AZ_CLEAR, data_org[AZ_name], label, (int*)0);

        /* If we're about to exit due to getting k==AZ_NEW_ADDRESS, let's try
           one other thing. See if the required memory can be found using
           type==scaling->mat_name. If not, then go ahead and crash. ABW
        */
        if (scaling != 0) {
          sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                               sizeof(double), AZ_ALLOC,
                                               scaling->mat_name, label, &k);
        }
      }

    if (k == AZ_NEW_ADDRESS) {
      (void) AZ_printf_err( "%sERROR: Previous scaling not found for matrix "
                     "with\ndata_org[AZ_name] = %d. Check\n"
                     "options[AZ_pre_calc]\n", yo,
                     data_org[AZ_name]);
      exit(-1);
    }
  }

  if ( (options[AZ_pre_calc] <= AZ_recalc) && (action == AZ_SCALE_MAT_RHS_SOL)) {
    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
      for (irow = 0; irow < N; irow++) {

        /* get row sum */

        j_last  = bindx[irow+1] - bindx[irow];
        bindx_row = bindx[irow];
        row_sum = fabs(val[irow]);

        for(jcol = 0; jcol < j_last; jcol++) {
          k        = bindx_row + jcol;
          row_sum += fabs(val[k]);
        }

        row_sum = row_sum * AZ_SGN(val[irow]);

        if (fabs(row_sum) < DBL_MIN) {
          (void) AZ_printf_err( "%sERROR: Row %d is all zero's\n", yo, irow);
          exit(-1);
        }

        sc_vec[irow] = 1.0 / row_sum;

        /* scale matrix row */

        val[irow] *= sc_vec[irow];
        for (jcol = 0; jcol < j_last; jcol++) {
          k       = bindx_row + jcol;
          val[k] *= sc_vec[irow];
        }
      }
    }

    else {
      m = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

      /* loop over the block rows */

      for (iblk_row = 0; iblk_row < m; iblk_row++) {

        /* find out how many rows are in this block row */

        num_blk_rows = rpntr[iblk_row+1] - rpntr[iblk_row];

        /* find out how many block columns are in this block row */

        num_col_blks = bpntr[iblk_row+1] - bpntr[iblk_row];

        /* loop over all the rows in this block row */

        for (irow = 0; irow < num_blk_rows; irow++) {
          I = rpntr[iblk_row] + irow;     /* true matrix row number */

          /* loop over all the column blocks in this block row */

          for (jblk_col = 0; jblk_col < num_col_blks; jblk_col++) {

            /* find out which column block this is */

            icol_blk = bindx[iblk];
            indx_ptr = indx[iblk++];

            /* find out how many columns are in this block */

            num_blk_cols = cpntr[icol_blk+1] - cpntr[icol_blk];

            /* loop over all the columns in this block */

            for (jcol = 0; jcol < num_blk_cols; jcol++) {
              J     = cpntr[icol_blk] + jcol;     /* true matrix column */
              t_index = indx_ptr + jcol*num_blk_rows + irow;

              /* diagonal entry => get sign */

              if (I == J) sign = AZ_SGN(val[t_index]);

              row_sum += fabs(val[t_index]);
            }
          }

          /* reset the block counter */

          iblk -= num_col_blks;

          if ( fabs(sign) < (1.0 - sqrt(DBL_EPSILON)) ) {
            (void) AZ_printf_err( "%sERROR: sign not set => no diagonal "
                           "entry.\n         inside row_sum.\n", yo);
            exit(-1);
          }

          if (fabs(row_sum) == 0.0) {
            (void) AZ_printf_err("%sERROR: row %d is all 0's.\n", yo, I);
            exit(-1);
          }

          inv_row_sum = sign / row_sum;
          sc_vec[I]   = inv_row_sum;
          row_sum     = sign = 0.0;

          /* scale the matrix */

          for (jblk_col = 0; jblk_col < num_col_blks; jblk_col++) {

            /* find out which column block this is */

            icol_blk = bindx[iblk];
            indx_ptr = indx[iblk++];

            /* find out how many columns are in this block */

            num_blk_cols = cpntr[icol_blk+1] - cpntr[icol_blk];

            /* loop over all the columns in this block */

            for (jcol = 0; jcol < num_blk_cols; jcol++) {
              J           = cpntr[icol_blk] + jcol;     /* true matrix column */
              t_index       = indx_ptr + jcol * num_blk_rows + irow;
              val[t_index] *= inv_row_sum;
            }
          }

          /* reset the block counter */

          iblk -= num_col_blks;
        }

        /* last row in this block row => offset the correction above */

        iblk += num_col_blks;
      }
    }
  }
  if ( (action == AZ_SCALE_MAT_RHS_SOL) || (action == AZ_SCALE_RHS) ) {
       for (irow = 0; irow < N; irow++) b[irow] *= sc_vec[irow];
  }
  if ( action == AZ_INVSCALE_RHS)  {
       for (irow = 0; irow < N; irow++) b[irow] /= sc_vec[irow];
  }

} /* AZ_row_sum_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sym_diagonal_scaling(int action, AZ_MATRIX *Amat, 
		        double b[], double x[], 
                        int options[], 
			int proc_config[], struct AZ_SCALING *scaling)

/*******************************************************************************

  Routine to symmetrically diagonally scale sparse matrix problem; Note: this
  scales the entire matrix problem Ax = b, the routine sym_rescale must be used
  to transform solution back to recover solution to original problem.

  Author:          John N. Shadid, SNL, 1421 (MSR format)
  =======          Lydie Prevost, SNL 1422 (VBR format) 

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                    Aztec Users Guide).

  bindx:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file  Aztec Users Guide).

  b:               Right hand side of linear system.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see  Aztec Users Guide).

  options:         Determines specific solution method and other parameters.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  x:               Current solution vector.

*******************************************************************************/

{

  /* local variables */

  register int j, k, irow, icol;
  int          N, m;
  int          j_last, bindx_row, i;
  double       *sc_vec;
  int count;
  char         label[80];

  char        *yo = "AZ_sym_diagonal_scaling: ";
  int         *indx, *bindx, *rpntr, *cpntr, *bpntr, *data_org;
  double      *val;


  /**************************** execution begins ******************************/

  val  = Amat->val;
  indx = Amat->indx;
  bindx = Amat->bindx;
  rpntr = Amat->rpntr;
  cpntr = Amat->cpntr;
  bpntr = Amat->bpntr;
  data_org = Amat->data_org;
  if (action == AZ_INVSCALE_SOL) {
     AZ_sym_reinvscale_sl(x, data_org, options, proc_config, scaling); return;
  }
  if (action == AZ_SCALE_SOL) {
     AZ_sym_rescale_sl(x, data_org, options, proc_config, scaling); return;
  }

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];
  m    = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];



  sprintf(label,"sc_vec%d",options[AZ_recursion_level]);
  sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                       sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], label, &i);
  scaling->action = AZ_left_and_right_scaling;
  if (((action == AZ_SCALE_RHS) || (action == AZ_INVSCALE_RHS) || 
       (options[AZ_pre_calc] >= AZ_reuse)) && (i == AZ_NEW_ADDRESS)) {

    if (i == AZ_NEW_ADDRESS) {
      AZ_manage_memory(0, AZ_CLEAR, data_org[AZ_name], label, (int*)0);

      /* If we're about to exit due to getting i==AZ_NEW_ADDRESS, let's try
         one other thing. See if the required memory can be found using
         type==scaling->mat_name. If not, then go ahead and crash. ABW
      */
      if (scaling != 0) {
        sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                             sizeof(double), AZ_ALLOC,
                                             scaling->mat_name, label, &i);
      }
    }

    if (i == AZ_NEW_ADDRESS) {
      (void) AZ_printf_err( "%sERROR: Previous scaling not found for matrix "
                   "with\ndata_org[AZ_name] = %d. Check\n"
                   "options[AZ_pre_calc]\n\n", yo,
                   data_org[AZ_name]);
      exit(-1);
    }
  }



  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {

    /***** MSR Matrix *****/


    if ( (options[AZ_pre_calc] <= AZ_recalc) && (action == AZ_SCALE_MAT_RHS_SOL)) {
      for (irow = 0; irow < N; irow++) {

        /* scale matrix */

        j_last  = bindx[irow+1] - bindx[irow];
        bindx_row = bindx[irow];

        if (fabs(val[irow]) < DBL_MIN) {
          (void) AZ_printf_err( "%sERROR: diagonal of row %d is zero\n", yo,
                         irow);
          exit(-1);
        }

        sc_vec[irow] = 1.0 / sqrt(fabs(val[irow]));

        for (j = 0; j < j_last; j++) {
          k       = bindx_row + j;
          val[k] *= sc_vec[irow];
        }
        val[irow] *= sc_vec[irow];
      }

      /* do right diagonal scaling */

      AZ_exchange_bdry(sc_vec, data_org, proc_config);

      /* index through rows of matrix */

      for (irow = 0; irow < N; irow++) {
        val[irow] *= sc_vec[irow];

        j_last     = bindx[irow+1] - bindx[irow];
        bindx_row    = bindx[irow];

        for (j = 0; j < j_last; j++) {
          k       = bindx_row + j;
          val[k] *= sc_vec[bindx[k]];
        }
      }
    }

  }

  else {

    /***** VBR Matrix *****/
   
    /* find the diagonal terms */

    if ( (options[AZ_pre_calc] <= AZ_recalc) && (action == AZ_SCALE_MAT_RHS_SOL)) {

     /* index through block rows of matrix */
   
         for (irow = 0; irow < m; irow++) {
   
         /* for all the blocks in a row */
   
         for (k = bpntr[irow]; k < bpntr[irow+1]; k++) {
            icol = bindx[k];
   
            count = 0;
   
            for (j = cpntr[icol]; j < cpntr[icol+1]; j++) {
               for (i = rpntr[irow]; i < rpntr[irow+1]; i++) {
                  if ( icol == irow && i ==j ){
                     sc_vec[i] =  1.0 / sqrt(fabs(val[indx[k]+count]));
                  }
                  count++;
               }
            }
   
         }
      }


     /* do left and right diagonal scaling */
   
      AZ_exchange_bdry(sc_vec, data_org, proc_config);
   
     /* index through rows of matrix */
   
      for (irow = 0; irow < m; irow++) {
   
         /* for all the blocks in a row */
   
         for (k = bpntr[irow]; k < bpntr[irow+1]; k++) {
            icol = bindx[k];
   
            count = 0;
   
            for (j = cpntr[icol]; j < cpntr[icol+1]; j++) {
               for (i = rpntr[irow]; i < rpntr[irow+1]; i++) {
                  val[indx[k]+count] *=  sc_vec[i] * sc_vec[j] ;
                  count++;
               }
            }
   
         }
      }
   
    }
  }

  /* rescale right hand side and solution */

  if ( (action == AZ_SCALE_RHS) ) {
       for (irow = 0; irow < N; irow++) b[irow] *= sc_vec[irow];
  }
  if ( action == AZ_INVSCALE_RHS)  {
       for (irow = 0; irow < N; irow++) b[irow] /= sc_vec[irow];
  }
  if (action == AZ_SCALE_MAT_RHS_SOL) {
       for (irow = 0; irow < N; irow++) b[irow] *= sc_vec[irow];
       for (irow = 0; irow < N; irow++) x[irow] /= sc_vec[irow];
  }

} /* sym_diagonal_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sym_row_sum_scaling(int action, AZ_MATRIX *Amat, 
                              double b[],
                              double x[], int options[], 
			      int proc_config[], struct AZ_SCALING *scaling)

/*******************************************************************************

  Routine to symmetrically diagonally scale sparse matrix problem; Note: this
  scales the entire matrix problem Ax = b, the routine sym_rescale must be used
  to transform solution back to recover solution to original problem.

  Author:          John N. Shadid, SNL, 1421 (MSR format)
  =======          David Day,      SNL, 9214 (VBR format)


  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                    Aztec Users Guide).

  bindx:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file  Aztec Users Guide).

  b:               Right hand side of linear system.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see  Aztec Users Guide).

  options:         Determines specific solution method and other parameters.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  x:               Current solution vector.

  Implementation note:  Symmetric row scaling initially matched the sign of
  the diagonal entries, and terminated if a diagonal vanished.  This feature
  (i.e. mistake ) has been discontinued, and zero rows are scaled by one.
*******************************************************************************/
{

  /* local variables */

  register int indx_ptr, row, col, blkrow, blkcol, locblk, iblk = 0;
  int          i, j, j_last, bindx_row, k;
  int          rblksz, numblk, cblksz, index;
  int          N;   /* number of unknowns  (rows) */
  int          m;   /* number of blocks (c.f. vbr) */
  double       row_sum = 0.0, *sc_vec, one = 1.0;
  char        *yo = "AZ_sym_row_sum_scaling: ";
  char         label[80];
  int         *indx, *bindx, *rptr, *cptr, *bptr, *data_org;
  double      *val;

  /**************************** execution begins ******************************/

  val = Amat->val;
  indx = Amat->indx;
  bindx = Amat->bindx;
  rptr = Amat->rpntr;
  cptr = Amat->cpntr;
  bptr = Amat->bpntr;
  data_org = Amat->data_org;
  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  if (action == AZ_INVSCALE_SOL) {
     AZ_sym_reinvscale_sl(x, data_org, options, proc_config, scaling); return;
  }
  if (action == AZ_SCALE_SOL) {
     AZ_sym_rescale_sl(x, data_org, options, proc_config, scaling); return;
  }

  sprintf(label,"sc_vec%d",options[AZ_recursion_level]);
  sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                       sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], label, &i);
  scaling->action = AZ_left_and_right_scaling;

  if ((options[AZ_pre_calc] >= AZ_reuse) && (i == AZ_NEW_ADDRESS)) {

    if (i == AZ_NEW_ADDRESS) {
        AZ_manage_memory(0, AZ_CLEAR, data_org[AZ_name], label, (int*)0);

      /* If we're about to exit due to getting i==AZ_NEW_ADDRESS, let's try
         one other thing. See if the required memory can be found using
         type==scaling->mat_name. If not, then go ahead and crash. ABW
      */
      if (scaling != 0) {
        sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                             sizeof(double), AZ_ALLOC,
                                             scaling->mat_name, label, &i);
      }
    }

    if (i == AZ_NEW_ADDRESS) {
      (void) AZ_printf_err( "%sERROR: Previous scaling not found for matrix "
                     "with\ndata_org[AZ_name] = %d. Check\n"
                     "options[AZ_pre_calc]\n", yo,
                     data_org[AZ_name]);
      exit(-1);
    }
  }
  if ( (options[AZ_pre_calc] <= AZ_recalc) && (action == AZ_SCALE_MAT_RHS_SOL)) {
    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {

	/***** MSR Matrix *****/

    for(row = 0; row < N; row++) {

      /* get row sum */

      j_last  = bindx[row+1] - bindx[row];
      bindx_row = bindx[row];
      row_sum = fabs(val[row]);

      for (col = 0; col < j_last; col++) {
        k        = bindx_row + col;
        row_sum += fabs(val[k]);
      }

      if (fabs(row_sum) < DBL_MIN) row_sum = one;

      sc_vec[row] = one/ sqrt(fabs(row_sum));

      /* scale matrix row */

      val[row] *= sc_vec[row];

      for (col = 0; col < j_last; col++) {
        k       = bindx_row + col;
        val[k] *= sc_vec[row];
      }
    }

    AZ_exchange_bdry(sc_vec, data_org, proc_config);

    /* do right diagonal scaling */
    /* index through rows of matrix */

    for (row = 0; row < N; row++) {
      val[row] *= sc_vec[row];
      j_last     = bindx[row+1] - bindx[row];
      bindx_row    = bindx[row];

      for (j = 0; j < j_last; j++) {
        k       = bindx_row + j;
        val[k] *= sc_vec[bindx[k]];
      }
    }
  } else {
      
      /***** VBR Matrix *****/

      m = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

      for (blkrow = 0; blkrow < m; blkrow++) {
        rblksz = rptr[blkrow+1] - rptr[blkrow];
        numblk = bptr[blkrow+1] - bptr[blkrow];
        for (row = 0; row < rblksz; row++) {
          i = rptr[blkrow] + row;
/*        
          Debugging tool
          if( i >= N  ){
            AZ_printf_out("i %d > N %d   m %d > blkrow %d  and  %d > %d\n",
                    i, N, m, blkrow, rblksz, row);
          }
*/
          for (blkcol = 0; blkcol < numblk; blkcol++) {
            locblk = bindx[iblk];
            indx_ptr = indx[iblk++];
            cblksz = cptr[locblk+1] - cptr[locblk];
            for (col = 0; col < cblksz; col++) {
              index = indx_ptr + col*rblksz + row;
              row_sum += fabs(val[index]);
/*        
              Debugging tool
              if( index >= nnzers ){
                AZ_printf_out("N_nz %d < index %d\n", nnzers, index);
              }
*/




            }
          }
          if (fabs(row_sum) < DBL_MIN){ 
/*        
              Debugging tool
              AZ_printf_out("Zed: i= %d   %d > %d  and  %d > %d\n",
                      i, numblk, blkcol, rblksz, row);
*/
            row_sum = one;
          }
          sc_vec[i]   = one/ row_sum;
          row_sum     = 0.0;
          iblk -= numblk;
        }
        iblk += numblk;
      }
      AZ_exchange_bdry(sc_vec, data_org, proc_config);

      /* Symmetricly scale the matrix */

      iblk = 0;
      for (blkrow = 0; blkrow < m; blkrow++) {
        rblksz = rptr[blkrow+1] - rptr[blkrow];
        numblk = bptr[blkrow+1] - bptr[blkrow];
        for (row = 0; row < rblksz; row++) {
          i = rptr[blkrow] + row;  
          row_sum = sc_vec[i];     /* Single look up */
          for (blkcol = 0; blkcol < numblk; blkcol++) {
            locblk = bindx[iblk];
            indx_ptr = indx[iblk++];
            cblksz = cptr[locblk+1] - cptr[locblk];
            for (col = 0; col < cblksz; col++) {
              j = cptr[locblk] + col;
              index = indx_ptr + col*rblksz + row;
              val[index] *= (row_sum * sc_vec[j]);
            }
          }
          iblk -= numblk;
        }
        iblk += numblk;
      }
    } /* VBR */
  }

  if ( (action == AZ_SCALE_RHS) ) {
       for (row = 0; row < N; row++) b[row] *= sc_vec[row];
  }
  if ( action == AZ_INVSCALE_RHS)  {
       for (row = 0; row < N; row++) b[row] /= sc_vec[row];
  }
  if (action == AZ_SCALE_MAT_RHS_SOL) {
       for (row = 0; row < N; row++) b[row] *= sc_vec[row];
       for (row = 0; row < N; row++) x[row] /= sc_vec[row];
  }
} /* AZ_sym_row_sum_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_equil_scaling(int action, AZ_MATRIX *Amat,
                              double b[],
                              double x[], int options[],
                              int proc_config[], struct AZ_SCALING *scaling)

/*******************************************************************************

  Symmetricly diagonally scale a sparse matrix.
  Use sym_rescale to recover the solution.  The
  scaling that equilibriates the matrix results
  from accumulating a few iterations of symmetric
  row sum scaling.

  Author:          David Day,      SNL, 9214


  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                   the Aztec Users Guide).

  bindx:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   the Aztec Users Guide).

  b:               Right hand side of linear system.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see the Aztec Users Guide).

  options:         Determines specific solution method and other parameters.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  x:               Current solution vector.

*******************************************************************************/
{

  /* local variables */

  register int indx_ptr, row, col, blkrow, blkcol, locblk, iblk = 0;
  int          i, j, j_last, bindx_row;
  int          rblksz, numblk, cblksz, index;
  int          m, k, N, Npe, iteration;
  double      *sc_vec, *rowsum, *colsum, one = 1.0, zed= 0.0, scaledval;
  double       maxrow, minrow;
  double       maxcol, mincol;
  char        *yo = "AZ_equilibration: ";
  char         label[80];
  int         *indx, *bindx, *rptr, *cptr, *bptr, *data_org;
  int          p;
  double      *val;

  /******************************** ExEcutE ***********************************/

  val = Amat->val;
  indx = Amat->indx;
  bindx = Amat->bindx;
  rptr = Amat->rpntr;
  cptr = Amat->cpntr;
  bptr = Amat->bpntr;
  data_org = Amat->data_org;
  N        = data_org[AZ_N_internal] + data_org[AZ_N_border];
  Npe      = N + data_org[AZ_N_external];
  p          = proc_config[AZ_node];

  if (action == AZ_INVSCALE_SOL) {
     AZ_sym_reinvscale_sl(x, data_org, options, proc_config, scaling); return;
  }
  if (action == AZ_SCALE_SOL) {
     AZ_sym_rescale_sl(x, data_org, options, proc_config, scaling); return;
  }

  sprintf(label,"sc_vec%d",options[AZ_recursion_level]);
  sc_vec = (double *) AZ_manage_memory( Npe * sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], label, &i);
  scaling->action = AZ_left_and_right_scaling;

  if ((options[AZ_pre_calc] >= AZ_reuse) && (i == AZ_NEW_ADDRESS)) {

    if (i == AZ_NEW_ADDRESS) {
      AZ_manage_memory(0, AZ_CLEAR, data_org[AZ_name], label, (int*)0);

      /* If we're about to exit due to getting i==AZ_NEW_ADDRESS, let's try
         one other thing. See if the required memory can be found using
         type==scaling->mat_name. If not, then go ahead and crash. ABW
      */
      if (scaling != 0) {
        sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                             sizeof(double), AZ_ALLOC,
                                             scaling->mat_name, label, &i);
      }
    }

    if (i == AZ_NEW_ADDRESS) {
      (void) AZ_printf_err( "%sERROR: Previous scaling not found for matrix "
                     "with\ndata_org[AZ_name] = %d. Check\n"
                     "options[AZ_pre_calc]\n", yo, data_org[AZ_name]);
      exit(-1);
    }
  }
  if ( (options[AZ_pre_calc] <= AZ_recalc) && (action == AZ_SCALE_MAT_RHS_SOL)) {

    rowsum = (double *) malloc(Npe*sizeof(double));
    colsum = (double *) malloc(Npe*sizeof(double));
    if (colsum == NULL){
      perror("Error: Not enough space to create matrix");
    }
    else if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {

      if( p == 0 ){
         printf("Untested code.  Comments to dmday@sandia.gov\n");
      }
      /***** MSR Matrix *****/

      for(i= 0; i< Npe; i++) {
        rowsum[i] = zed;
        sc_vec[i] = one;
      }

      for(iteration=0; iteration<5; iteration++){
        for(i= 0; i< Npe; i++) {
          colsum[i] = zed;
        }
        scaledval = fabs(val[0])*sc_vec[0]*sc_vec[0];
        maxrow = zed;
        minrow = 1.e+40;

        for(row = 0; row < N; row++) {
          bindx_row = bindx[row];
          j_last  = bindx[row+1] - bindx_row;
          scaledval = fabs(val[row])*sc_vec[row]*sc_vec[row];

          /* Initialize sums to diagonal values */
          rowsum[row] = scaledval;
          colsum[row] = scaledval;

          for (j = 0; j < j_last; j++) {
            k        = bindx_row + j;
            col      = bindx[k];
            scaledval = fabs(val[k])*sc_vec[row]*sc_vec[col];
            rowsum[row] += scaledval;
            colsum[col] += scaledval;
          }

          if( rowsum[row] < minrow ) minrow = rowsum[row];
          if( rowsum[row] > maxrow ) maxrow = rowsum[row];
        }
        AZ_sum_bdry(colsum, data_org, proc_config);

        /* Update Scaling */
        for(row = 0; row < N; row++) {
          scaledval = rowsum[row] + colsum[row];
          if( scaledval > DBL_MIN){
            if (rowsum[row] < DBL_MIN){
              scaledval = colsum[row];
            } else if (colsum[row] < DBL_MIN){
              scaledval = rowsum[row];
            } else {
              scaledval = sqrt( rowsum[row] * colsum[row] );
            }
            sc_vec[row] /= sqrt(scaledval);
          }
        }
        AZ_exchange_bdry(sc_vec, data_org, proc_config);
      } /* iteration */

      /* Scale Matrix */
      for(row = 0; row < N; row++) {
        scaledval   = sc_vec[row]*sc_vec[row];
        val[row]   *= scaledval;
        bindx_row   = bindx[row];
        j_last      = bindx[row+1] - bindx_row;
        for (j = 0; j < j_last; j++) {
          k        = bindx_row + j;
          col      = bindx[k];
          scaledval = sc_vec[row]*sc_vec[col];
          val[k]   *= scaledval;
        }
      }
    } else {

      /***** VBR Matrix *****/

      for(i= 0; i< Npe; i++) {
        sc_vec[i] = one;
      }
      m = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

      for(iteration=0; iteration<5; iteration++){
        for(i= 0; i< Npe; i++) {
          rowsum[i] = zed;
          colsum[i] = zed;
        }
        scaledval = fabs(val[indx[0]])*sc_vec[0]*sc_vec[0];
        maxrow = zed;
        minrow = 1.e+40;

        iblk = 0;
        for (blkrow = 0; blkrow < m; blkrow++) {
          rblksz = rptr[blkrow+1] - rptr[blkrow];
          numblk = bptr[blkrow+1] - bptr[blkrow];
          for (row = 0; row < rblksz; row++) {
            i = rptr[blkrow] + row;
            for (blkcol = 0; blkcol < numblk; blkcol++) {
              locblk = bindx[iblk];
              indx_ptr = indx[iblk++];
              cblksz = cptr[locblk+1] - cptr[locblk];
              for (col = 0; col < cblksz; col++) {
                j = cptr[locblk] + col;
                index = indx_ptr + col*rblksz + row;
                scaledval = fabs(val[index])*sc_vec[i]*sc_vec[j];
                rowsum[i] += scaledval;
                colsum[j] += scaledval;
              }
            }
            iblk -= numblk;
            if( rowsum[i] < minrow ) minrow = rowsum[i];
            if( rowsum[i] > maxrow ) maxrow = rowsum[i];
          }
          iblk += numblk;
        }

        /* Complete column sums */
        AZ_sum_bdry(colsum, data_org, proc_config);

        /* Diagnostics */
        maxcol = colsum[0]; 
        mincol = colsum[0]; 
        for (col = 1; col < N; col++) {
          if( colsum[col] < mincol ) mincol = colsum[col];
          if( colsum[col] > maxcol ) maxcol = colsum[col];
        }
        if( p == 0 )
        printf("(%d) AZ_equil: enter sum %d    %e row %e    %e col %e\n",
                p,iteration, minrow, maxrow, mincol, maxcol);

        /* Update Scaling */
        for (blkrow = 0; blkrow < m; blkrow++) {
          rblksz = rptr[blkrow+1] - rptr[blkrow];
          numblk = bptr[blkrow+1] - bptr[blkrow];
          for (row = 0; row < rblksz; row++) {
            i = rptr[blkrow] + row;
            scaledval = rowsum[i] + colsum[i];
            if( scaledval > DBL_MIN){
              if (rowsum[i] < DBL_MIN){
                scaledval = colsum[i];
              } else if (colsum[i] < DBL_MIN){
                scaledval = rowsum[i];
              } else {
                scaledval = sqrt( rowsum[i] * colsum[i] );
              }
              sc_vec[i] /= sqrt(scaledval);
            }
          }
        }
        AZ_exchange_bdry(sc_vec, data_org, proc_config);
      } /* iteration */

      /* Symmetricly scale the matrix */
      iblk = 0;
      for (blkrow = 0; blkrow < m; blkrow++) {
        rblksz = rptr[blkrow+1] - rptr[blkrow];
        numblk = bptr[blkrow+1] - bptr[blkrow];
        for (row = 0; row < rblksz; row++) {
          i = rptr[blkrow] + row;
          scaledval = sc_vec[i];     /* Single look up */
          for (blkcol = 0; blkcol < numblk; blkcol++) {
            locblk = bindx[iblk];
            indx_ptr = indx[iblk++];
            cblksz = cptr[locblk+1] - cptr[locblk];
            for (col = 0; col < cblksz; col++) {
              j = cptr[locblk] + col;
              index = indx_ptr + col*rblksz + row;
              val[index] *= (scaledval * sc_vec[j]);
            }
          }
          iblk -= numblk;
        }
        iblk += numblk;
      }
    } /* VBR */
    free((void *) rowsum);
    free((void *) colsum);
  }

  if ( (action == AZ_SCALE_RHS) ) {
       for (row = 0; row < N; row++) b[row] *= sc_vec[row];
  }
  if ( action == AZ_INVSCALE_RHS)  {
       for (row = 0; row < N; row++) b[row] /= sc_vec[row];
  }
  if (action == AZ_SCALE_MAT_RHS_SOL) {
       for (row = 0; row < N; row++) b[row] *= sc_vec[row];
       for (row = 0; row < N; row++) x[row] /= sc_vec[row];
  }
} /* AZ_equil_scaling */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sym_rescale_sl(double x[], int data_org[], int options[],
	int proc_config[], struct AZ_SCALING *scaling)

/*******************************************************************************

  Routine to symmetrically diagonally rescale the sparse matrix problem;
  Note: this rescales the entire matrix problem Ax = b.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               Current solution vector.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see  Aztec Users Guide).

*******************************************************************************/

{

  /* local variables */

  register int i;
  int          N, k;
  double      *sc_vec;
  char        *yo = "AZ_sym_rescale_sl: ";
  char        label[80];

  /**************************** execution begins ******************************/

  if ((data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
      (data_org[AZ_matrix_type] != AZ_VBR_MATRIX) ) {
     (void) AZ_printf_err("%sWARNING: Matrix type is neither MSR nor VBR\n",
                    yo);
     return;
  }



  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  sprintf(label,"sc_vec%d",options[AZ_recursion_level]);
  sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                       sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], label, &k);
  scaling->action = AZ_left_and_right_scaling;

  if (k == AZ_NEW_ADDRESS) {
    (void) AZ_printf_err( "%sWARNING: Scaling vector not found: "
                   "Not rescaling solution\n", yo);
    return;
  }

  for (i = 0; i < N; i++) x[i] = x[i] / sc_vec[i];

  AZ_exchange_bdry(x, data_org, proc_config);

} /* AZ_sym_rescale_sl */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sym_reinvscale_sl(double x[], int data_org[], int options[],
	int proc_config[], struct AZ_SCALING *scaling)

/*******************************************************************************

  Routine to symmetrically diagonally rescale the sparse matrix problem;
  Note: this rescales the entire matrix problem Ax = b.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               Current solution vector.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see  Aztec Users Guide).

*******************************************************************************/

{

  /* local variables */

  register int i;
  int          N, k;
  double      *sc_vec;
  char        *yo = "AZ_sym_rescale_sl: ";
  char        label[80];

  /**************************** execution begins ******************************/

  if ((data_org[AZ_matrix_type] != AZ_MSR_MATRIX) &&
      (data_org[AZ_matrix_type] != AZ_VBR_MATRIX) ) {
     (void) AZ_printf_err("%sWARNING: Matrix type is neither MSR nor VBR\n",
                    yo);
     return;
  }



  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  sprintf(label,"sc_vec%d",options[AZ_recursion_level]);
  sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                       sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], label, &k);
  scaling->action = AZ_left_and_right_scaling;

  if (k == AZ_NEW_ADDRESS) {
    (void) AZ_printf_err( "%sWARNING: Scaling vector not found: "
                   "Not rescaling solution\n", yo);
    return;
  }

  for (i = 0; i < N; i++) x[i] = x[i] * sc_vec[i];

  AZ_exchange_bdry(x, data_org, proc_config);

} /* AZ_sym_rescale_sl */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void calc_blk_diag_Chol(double *val, int *indx, int *bindx, int *rpntr,
                               int *cpntr, int *bpntr, double *L, int *d_indx,
                               int *d_bindx, int *d_rpntr, int *d_bpntr,
                               int *data_org)

/*******************************************************************************

  Routine to calculate the Cholesky factors of the block-diagonal portion of the
  sparse matrix in 'val' and the associated integer pointer vectors. This is
  used for scaling and/or preconditioning.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                    Aztec Users Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file  Aztec Users Guide).

  L:               Vector containing the upper Cholesky factors of the diagonal
                   blocks.

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
                   (see  Aztec Users Guide).

*******************************************************************************/

{

  /* local variables */

  register int i, j, iblk_row, jblk, icount = 0, iblk_count = 0, ival;
  int          m1, n1, itemp, iL;
  int          m;
  int          bpoff, idoff;
  int          info;
  char        *yo = "calc_blk_diag_Chol: ", *uplo = "L";

  /**************************** execution begins ******************************/

  m = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

  if (m == 0) return;

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
          (void) AZ_printf_err( "%sERROR: diagonal blocks are not square.\n",
                         yo);
          exit(-1);
        }
        else {

          /* fill the vectors */

          d_indx[iblk_count]  = icount;
          d_rpntr[iblk_count] = rpntr[iblk_row];
          d_bpntr[iblk_count] = iblk_row;
          d_bindx[iblk_count] = iblk_row;

          for (i = 0; i < itemp; i++) L[icount++] = val[ival + i];

          /* Compute the Cholesky factors for this block */

          iL = d_indx[d_bpntr[iblk_row] - *d_bpntr] - *d_indx;
          DPOTRF_F77(CHAR_MACRO(uplo[0]), &m1, L+iL, &m1, &info);

          if (info < 0) {
            (void) AZ_printf_err( "%sERROR: argument %d is illegal\n", yo,
                           -info);
            exit(-1);
          }
          else if (info > 0) {
            (void) AZ_printf_err( "%sERROR: the factorization has produced a "
                           "singular L with L[%d][%d] being exactly zero\n",
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

} /* calc_blk_diag_Chol */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

#ifdef next_version
void AZ_sym_rescale_vbr(double x[], int data_org[], int options[])

/*******************************************************************************

  Routine to symmetrically block-diagonally rescale the sparse matrix problem;
  Note: this rescales the entire matrix problem Ax = b.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               Current solution vector.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see  Aztec User)'s Guide.

*******************************************************************************/

{

  /* local variables */

  int          N, k;
  double      *sc_vec;
  char        *yo = "AZ_sym_rescale_vbr: ";

  /**************************** execution begins ******************************/

  if (data_org[AZ_matrix_type] != AZ_VBR_MATRIX) return;

  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  sprintf(label,"sc_vec%d",options[AZ_recursion_level]);
  sc_vec = (double *) AZ_manage_memory((N + data_org[AZ_N_external]) *
                                       sizeof(double), AZ_ALLOC,
                                       data_org[AZ_name], label, &k);
  scaling->action = AZ_left_and_right_scaling;

  if (k == AZ_NEW_ADDRESS) {
    (void) AZ_printf_err( "%sWARNING: Scaling vector not found - "
                   "not rescaling solution\n", yo);
    return;
  }
  /*
    for (i = 0; i < N; i++) x[i] = x[i] * sc_vec[i];

    AZ_exchange_bdry(x, data_org, proc_config);
    */
} /* AZ_sym_rescale_vbr */
#endif

struct AZ_SCALING *AZ_scaling_create()
{
   struct AZ_SCALING *temp;

   temp = (struct AZ_SCALING *) AZ_allocate(sizeof(struct AZ_SCALING));
   if (temp == NULL) {
      AZ_printf_err("AZ_scaling_create: Not enough space\n");
      exit(1);
   }
   temp->A_norm = 0.0;
   temp->action = AZ_none;
   temp->mat_name = 0;
   temp->scaling_opt = AZ_none;
   temp->scale = 0;
   temp->scaling_data = 0;
   return temp;
} /* AZ_scaling_create() */

void AZ_scaling_destroy(struct AZ_SCALING **temp)
{
   if (*temp != NULL) AZ_free(*temp);
   *temp = NULL;
} /* AZ_scaling_destroy */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void AZ_sum_bdry(double x[], int data_org[], int proc_config[])

/*******************************************************************************

  Sum the external values to the corresponding border value. AZ_sum_bdry is 
  derived from the Schwarz preconditioning functions developed by Ray Tuminaro.

  Author:          David Day, SNL
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters

  x:               Vector of unknowns defined on the current processor.

*******************************************************************************/

{
  static int type = 0; /* avoid AZ_sys_msg_type bug */
  int name, total, i, j, count=0, st, from, N_unpadded;
  int nwords;
  MPI_AZRequest request[AZ_MAX_NEIGHBORS];  /* Message handle */
  double *buffer;                    /* for incoming messages */
  N_unpadded = data_org[AZ_N_internal] + data_org[AZ_N_border];
  name       = data_org[AZ_name];
  ++type;
  type       = type % AZ_NUM_MSGS;

  /* Send the external points to neighbors */
  total = 0;
  for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ )
     total += data_org[AZ_send_length+i];
  /*
  AZ_printf_out("(%d) mtype %d n_mes %d total words %d\n",p, type, data_org[AZ_N_neigh],total);
   */
  buffer= (double *) AZ_manage_memory(total*sizeof(double), AZ_ALLOC,
                     name, "temp in combine", &i);
  for ( i = 0 ; i < total ; i++ ) buffer[i] = 0.;

  /* Post receives */
  count = 0;
  for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ ) {
     from = data_org[AZ_neighbors+i];
     (void) mdwrap_iread((void *) &(buffer[count]),
                  sizeof(double)*data_org[AZ_send_length+i],
                  &from, &type, request+i);
     count += data_org[AZ_send_length+i];
     /*
      * AZ_printf_out("(%d) post receive %d from %d count %d\n",p, i, from, count);
      */
  }

  /* Send messages */
  count = N_unpadded;
  for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ ) {
     (void) mdwrap_write((void *) &(x[count]), data_org[AZ_rec_length+i]*
                     sizeof(double), data_org[AZ_neighbors+i], type, &st);
     count += data_org[AZ_rec_length+i];
  }
 
  /* Receive messages and sum recvd values to the send list */
  count = 0;
  for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ ) {
     from = data_org[AZ_neighbors+i];
     nwords = data_org[AZ_send_length+i];
     (void) mdwrap_wait((void *) &(buffer[count]),
                  sizeof(double)*nwords, &from, &type, &st,request+i);
     count += nwords;
     /*
      * AZ_printf_out("(%d) receive %d from %d count %d\n",p, i, from, count);
      */
  }
  for ( j = 0 ; j < total; j++ )
       x[ data_org[AZ_send_list + j] ] += buffer[j];
} /* AZ_sum_bdry */
