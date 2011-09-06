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
#include <math.h>
#include <stdio.h>
#include "az_aztec.h"
#include "az_blas_wrappers.h"

extern int az_iterate_id;

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_polynomial_expansion( double z[], int options[], int proc_config[], 
                              AZ_PRECOND *precond )
/*******************************************************************************

  Uses a Neuman series expansion to approximate the inverse of a matrix. The
  series expansion is in terms of (I - A/omega) where I is the identity, A the
  matrix for which the inverse is being approximated, and omega is a scaling
  factor (omega >= || A || / 2 , Wong and Jiang (1989) or the diagonal element
  if it is a constant). If power = 0 then diagonal scaling is performed. If
  power < 0 then an unparameterized expansion is used. If power > 0 then a
  parameterized expansion developed by a least squares method is used. This
  technique minimizes the L2 norm of the residual polynomial R(), on an evalue
  interval of [0,lambda_max] where lambda_max is an estimate of the largest
  evalue of A.(see Saad (1985)).

  This version assumes that diagonal scaling has been carried out on the entire
  set of equations.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  z:               On input, is the residual(rhs) of the set of equations.
                   On output is the result.

  options:         Determines specific solution method and other parameters.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  precond:         Structure used to represent the preocnditioner
                   (see az_aztec.h and Aztec User's Guide).

*******************************************************************************/

{

  /* local variables */

  int              param_flag, one = 1, j;
  register int     i, p;
  register double  cp;
  double           lambda_max;
  static double    c[15], inv_omega;
  int              N, power;
  double          *w, *poly_temp;
  int          *data_org, *bindx, *indx, *cpntr, *rpntr, *bpntr;
  double       *val;


  /**************************** execution begins ******************************/


  data_org = precond->Pmat->data_org;
  val =  precond->Pmat->val;
  bindx = precond->Pmat->bindx;
  cpntr = precond->Pmat->cpntr;
  indx = precond->Pmat->indx;
  rpntr = precond->Pmat->rpntr;
  bpntr = precond->Pmat->bpntr;
  N     = data_org[AZ_N_internal] + data_org[AZ_N_border];
  power = options[AZ_poly_ord];

  poly_temp = (double *) AZ_manage_memory(2*(N+data_org[AZ_N_external])*
                                          sizeof(double), AZ_ALLOC, AZ_SYS+az_iterate_id,
                                          "poly mem", &j);
  w         = &(poly_temp[N+data_org[AZ_N_external]]);

  if (options[AZ_precond] == AZ_Neumann ) param_flag = 0;
  else                                    param_flag = 1;

  if (options[AZ_pre_calc] < AZ_sys_reuse) {

    if (precond->Pmat->data_org[AZ_matrix_type] == AZ_USER_MATRIX) {
       lambda_max = precond->Pmat->matrix_norm;
       if (lambda_max < 0.0) {
           if (proc_config[AZ_node] == 0) {
               AZ_printf_err("Error: Matrix norm not given. Use ");
               AZ_printf_err("AZ_set_MATFREE_matrix_norm() to set it.\n");
           }
           exit(1);
       }
    }
    else if (precond->Pmat->data_org[AZ_matrix_type] == AZ_MSR_MATRIX ||
             precond->Pmat->data_org[AZ_matrix_type] == AZ_VBR_MATRIX ) {

       lambda_max = AZ_gmax_matrix_norm(val, indx, bindx, rpntr, cpntr, bpntr,
                                        proc_config, data_org);
   
       /* change sign of lambda_max if diagonal contains only negative values */
   
       AZ_change_sign(&lambda_max, val, indx, bindx, rpntr, cpntr, bpntr,
                      data_org);
   
    }
    inv_omega  = 1.0 / (0.55 * lambda_max);     /* 1.1*lambda_max/2 */
   
    if (param_flag)
      AZ_get_poly_coefficients(power, lambda_max, c, param_flag);
  }


  switch (param_flag) {
  case 0:                       /* Neumann series */
    DSCAL_F77(&N, &inv_omega, z, &one);
    DCOPY_F77(&N, z, &one, w, &one);

    for (p = power; p > 0; p--){
    precond->Pmat->matvec(z, poly_temp, precond->Pmat, proc_config);


      for (i = 0; i < N; i++)
        z[i] += w[i] - inv_omega * poly_temp[i];
    }
    break;

  case 1:                       /* least squares */

    /* initialization */

    DCOPY_F77(&N, z, &one, w, &one);
    DSCAL_F77(&N, c+power, z, &one);

    for (p = power - 1; p >= 0; p--) {
    precond->Pmat->matvec(z, poly_temp, precond->Pmat, proc_config);


      cp = *(c+p);
      for (i = 0; i < N; i++)
        z[i] = cp * w[i] + poly_temp[i];
    }
    break;

  default:
    if (proc_config[AZ_node] == 0) {
      (void) AZ_printf_err( "Error: invalid polynomial preconditioner\n"
                     "       options[AZ_precond] improperly set.\n");
    }
    exit(-1);
  }

} /* AZ_polynomial_expansion */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_get_poly_coefficients(int power, double b, double c[], int param_flag)

/*******************************************************************************

  Returns the coefficients for least-squares polynomial preconditioning.  From
  Saad, 1985.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  power:           Order of the polynomial.

  b:               Estimate for largest eigenvalue.

  c:               On output, requested coefficients.

  param_flag:      Determines if this is for Neumann series or least-squares
                   polynomial preconditioning.

*******************************************************************************/

{

  register int    i;
  register double temp;

  /* define temporary variable */

  temp = 4.0 / b;

  /* define least squares coefficients for polys ; evalue range [0,4] */
  /* Saad (1985) */

  if (param_flag) {

    switch(power) {
    case 0:
      c[0] = 1.0;
      break;
    case 1:
      c[0] = 5.0;      c[1] = -1.0;
      break;
    case 2:
      c[0] = 14.0;     c[1] = -7.0;     c[2] = 1.0;
      break;
    case 3:
      c[0] = 30.0;     c[1] = -27.0;    c[2] = 9.0; c[3] = -1.0;
      break;
    case 4:
      c[0] = 55.0;     c[1] = -77.0;    c[2] = 44.0;
      c[3] = -11.0;    c[4] = 1.0;
      break;
    case 5:
      c[0] = 91.0;     c[1] = -182.0;   c[2] = 156.0;
      c[3] = -65.0;    c[4] = 13.0;     c[5] = -1.0;
      break;
    case 6:
      c[0] = 140.0;    c[1] = -378.0;  c[2] = 450.0;
      c[3] = -275.0;   c[4] = 90.0;    c[5] = -15.0;
      c[6] = 1.0;
      break;
    case 7:
      c[0] = 204.0;    c[1] = -714.0;  c[2] = 1122.0;
      c[3] = -935.0;   c[4] = 442.0;   c[5] = -119.0;
      c[6] = 17.0;     c[7] = -1.0;
      break;
    case 8:
      c[0] = 285.0;    c[1] = -1254.0; c[2] = 2508.0;
      c[3] = -2717.0;  c[4] = 1729.0;  c[5] = -665.0;
      c[6] = 152.0;    c[7] = -19.0;   c[8] = 1.0;
      break;
    case 9:
      c[0] = 385.0;    c[1] = -2079.0; c[2] = 5148.0;
      c[3] = -7007.0;  c[4] = 5733.0;  c[5] = -2940.0;
      c[6] = 952.0;    c[7] = -189.0;  c[8] = 21.0;
      c[9] = -1.0;
      break;
    case 10:
      c[0] = 506.0;    c[1] = -3289.0; c[2] = 9867.0;
      c[3] = -16445.0; c[4] = 16744.0; c[5] = -10948.0;
      c[6] = 4692.0;   c[7] = -1311.0; c[8] = 230.0;
      c[9] = -23.0;    c[10] = 1.0;
      break;
    }

    /* transform coefficients to evalue interval [0,b] */

    for(i = 0; i <= power; i++)
      c[i] *= pow(temp, (double) i);
  }

  else
    for(i = 0; i <= power; i++)
      c[i] = 1.0;

} /* AZ_get_poly_coefficients */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
extern void AZ_funswill(int *);

void AZ_change_sign(double *lambda_max, double val[], int indx[], int bindx[],
                    int rpntr[], int cpntr[], int bpntr[], int data_org[])

/*******************************************************************************

  Figure out whether the diagonal has positive or negative numbers. If the
  diagonal contains only negative values, change the sign of lambda_max.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  lambda_max:      Estimate for the largest eigenvalue of the coefficient
                   matrix.

  val:             Array containing the nonzero entries of the matrix (see
                   Aztec User's Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

*******************************************************************************/

{

  /* local variables */

  int    pos, neg, i, kk, block_col, blk_size, j;
  double temp;

  /* external functions */


  /**************************** execution begins ******************************/

  pos = neg = 0;
  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
    for (i = 0; i < data_org[AZ_N_internal] + data_org[AZ_N_border]; i++) {
      if (val[i] > 0.0)      pos = 1;
      else if (val[i] < 0.0) neg = 1;
    }
  }

  else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
    for (i = 0; i < data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk]; i++) {
      for (kk = bpntr[i]; kk < bpntr[i+1]; kk++) {
        block_col = bindx[kk];

        if (block_col == i) {
          blk_size = cpntr[block_col+1] - cpntr[block_col];

          for (j = rpntr[i]; j < rpntr[i+1]; j++) {
            temp = val[indx[kk] + (blk_size+1) * (j-rpntr[i])];
            if (temp > 0.0)      pos = 1;
            else if (temp < 0.0) neg = 1;
          }
        }
      }

      AZ_funswill(&kk);         /* dummy so paragon compiler works right */
    }
  }

  if ((data_org[AZ_matrix_type] == AZ_MSR_MATRIX) ||
      (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)) {

    if ((pos == 0) && (neg == 0) && 
        (data_org[AZ_N_internal] + data_org[AZ_N_border] != 0) ) 
      (void) AZ_printf_err( "Warning: No nonzero matrix diagonal elements\n");

    if (pos + neg == 2) {
      (void) AZ_printf_err(
                     "Warning: Negative and positive matrix diagonal elements\n"
                     "         Better to use scaling with polynomial\n"
                     "         preconditioners in this case.\n");
    }
    else if (neg)
      *lambda_max = -(*lambda_max);
  }

} /* AZ_change_sign */
