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
#include <float.h>
#include "az_aztec.h"
#include "az_blas_wrappers.h"


extern void dvbr_sparax_basic(int m, double *val, int *bindx, int *rpntr,
                       int *cpntr, int *bpntr, double *b, double *c,
                       int exchange_flag, int *data_org, int *proc_config);

void dvbr_sparax_basic(int m, double *val, int *bindx, int *rpntr,
                       int *cpntr, int *bpntr, double *b, double *c,
                       int exchange_flag, int *data_org, int *proc_config)

/******************************************************************************

  c = Ab:
  Sparse (square) matrix-vector multiply, using the variable block row (VBR)
  data structure (A = val).

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  m:               Number of (block) rows in A.

  val:             Array containing the entries of the matrix. The matrix is
                   stored block-row-by-block-row. Each block entry is dense and
                   stored by columns (VBR).

  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   Aztec User's Guide).

  b:               Right hand side of linear system.

  c:               On output contains the solution to the linear system.

  exchange_flag:   Flag which controls call to AZ_exchange_bdry() (ignored in
                   serial implementation).

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

******************************************************************************/

{

  /* local variables */

  register double *x;
  register double *c_pntr;
  register int     iblk_row, j, jblk, iblk_size;
  int              m1, ib1, n1;
  int              bpoff, rpoff;
  int              ione = 1;
  int              irpntr, irpntr_next;
  int              ibpntr, ibpntr_next = 0;
  double           one = 1.0;
  double          *val_pntr;
  char            *N = "N";

  /**************************** execution begins *****************************/

  /* exchange boundary info */

  if (exchange_flag) AZ_exchange_bdry(b, data_org, proc_config);

  /* offset of the first block */

  bpoff = *bpntr;
  rpoff = *rpntr;

  /* zero the result vector */

  for (j = 0; j < rpntr[m] - rpoff; c[j++] = 0.0);

  val_pntr    = val;
  irpntr_next = *rpntr++;
  bpntr++;
  c          -= rpoff;

  /* loop over block rows */

  for (iblk_row = 0; iblk_row < m; iblk_row++) {

    irpntr      = irpntr_next;
    irpntr_next = *rpntr++;

    ibpntr      = ibpntr_next;
    ibpntr_next = *bpntr++ - bpoff;

    /* set result pointer */

    c_pntr      = c + irpntr;

    /* number of rows in the current row block */

    m1          = irpntr_next - irpntr;

    /* loop over each block in the current row block */

    for (j = ibpntr; j < ibpntr_next; j++) {
      jblk = *(bindx+j);

      /* the starting point column index of the current block */

      ib1 = *(cpntr+jblk);

      /* number of columns in the current block */

      n1 = cpntr[jblk + 1] - ib1;
      iblk_size = m1*n1;

      /****************** Dense matrix-vector multiplication *****************/

      /*
       * Get base addresses
       */

      x = b + ib1;


      /*
       * Special case the m1 = n1 = 1 case
       */

      if (iblk_size == 1)
        *c_pntr += *val_pntr * *x;

      else if (m1 == n1) {

        /*
         * Inline small amounts of work
         */

        switch (m1) {

        case 2:
          c_pntr[0] += val_pntr[0]*x[0] + val_pntr[2]*x[1];
          c_pntr[1] += val_pntr[1]*x[0] + val_pntr[3]*x[1];
          break;

        case 3:
          c_pntr[0] += val_pntr[0]*x[0] + val_pntr[3]*x[1] + val_pntr[6]*x[2];
          c_pntr[1] += val_pntr[1]*x[0] + val_pntr[4]*x[1] + val_pntr[7]*x[2];
          c_pntr[2] += val_pntr[2]*x[0] + val_pntr[5]*x[1] + val_pntr[8]*x[2];
          break;

        case 4:
          c_pntr[0] += val_pntr[0]*x[0] + val_pntr[4]*x[1] + val_pntr[8] *x[2]
            + val_pntr[12]*x[3];
          c_pntr[1] += val_pntr[1]*x[0] + val_pntr[5]*x[1] + val_pntr[9] *x[2]
            + val_pntr[13]*x[3];
          c_pntr[2] += val_pntr[2]*x[0] + val_pntr[6]*x[1] + val_pntr[10]*x[2]
            + val_pntr[14]*x[3];

          c_pntr[3] += val_pntr[3]*x[0] + val_pntr[7]*x[1] + val_pntr[11]*x[2]
            + val_pntr[15]*x[3];
          break;

        case 5:
          c_pntr[0] += val_pntr[0]*x[0] + val_pntr[5]*x[1] + val_pntr[10]*x[2]
            + val_pntr[15]*x[3] + val_pntr[20]*x[4];
          c_pntr[1] += val_pntr[1]*x[0] + val_pntr[6]*x[1] + val_pntr[11]*x[2]
            + val_pntr[16]*x[3] + val_pntr[21]*x[4];
          c_pntr[2] += val_pntr[2]*x[0] + val_pntr[7]*x[1] + val_pntr[12]*x[2]
            + val_pntr[17]*x[3] + val_pntr[22]*x[4];
          c_pntr[3] += val_pntr[3]*x[0] + val_pntr[8]*x[1] + val_pntr[13]*x[2]
            + val_pntr[18]*x[3] + val_pntr[23]*x[4];
          c_pntr[4] += val_pntr[4]*x[0] + val_pntr[9]*x[1] + val_pntr[14]*x[2]
            + val_pntr[19]*x[3] + val_pntr[24]*x[4];
          break;

        case 6:
          c_pntr[0] += val_pntr[0]*x[0] + val_pntr[6] *x[1] + val_pntr[12]*x[2]
            + val_pntr[18]*x[3] + val_pntr[24]*x[4] + val_pntr[30]*x[5];
          c_pntr[1] += val_pntr[1]*x[0] + val_pntr[7] *x[1] + val_pntr[13]*x[2]
            + val_pntr[19]*x[3] + val_pntr[25]*x[4] + val_pntr[31]*x[5];
          c_pntr[2] += val_pntr[2]*x[0] + val_pntr[8] *x[1] + val_pntr[14]*x[2]
            + val_pntr[20]*x[3] + val_pntr[26]*x[4] + val_pntr[32]*x[5];
          c_pntr[3] += val_pntr[3]*x[0] + val_pntr[9] *x[1] + val_pntr[15]*x[2]
            + val_pntr[21]*x[3] + val_pntr[27]*x[4] + val_pntr[33]*x[5];
          c_pntr[4] += val_pntr[4]*x[0] + val_pntr[10]*x[1] + val_pntr[16]*x[2]
            + val_pntr[22]*x[3] + val_pntr[28]*x[4] + val_pntr[34]*x[5];
          c_pntr[5] += val_pntr[5]*x[0] + val_pntr[11]*x[1] + val_pntr[17]*x[2]
            + val_pntr[23]*x[3] + val_pntr[29]*x[4] + val_pntr[35]*x[5];
          break;

        default:

          /*
           * For most computers, a really well-optimized assembly-coded level 2
           * blas for small blocks sizes doesn't exist.  It's better to
           * optimize your own version, and take out all the overhead from the
           * regular dgemv call.  For large block sizes, it's also a win to
           * check for a column of zeroes; this is what dgemv_ does.  The
           * routine dgemvnsqr_() is a fortran routine that contains optimized
           * code for the hp, created from the optimizing preprocessor. Every
           * workstation will probably have an entry here eventually, since
           * this is a key optimization location.
           */

/* #ifdef AZ_PA_RISC */
/*           dgemvnsqr_(&m1, val_pntr, x, c_pntr); */
/* #else */
          if (m1 < 10)
            AZ_dgemv2(m1, n1, val_pntr, x, c_pntr);
          else
            DGEMV_F77(CHAR_MACRO(N[0]), &m1, &n1, &one, val_pntr, &m1, x, &ione, &one, c_pntr,
                   &ione);
/* #endif */

        }
      }

      /* nonsquare cases */

      else {
        if (m1 < 10)
          AZ_dgemv2(m1, n1, val_pntr, x, c_pntr);
        else
          DGEMV_F77(CHAR_MACRO(N[0]), &m1, &n1, &one, val_pntr, &m1, x, &ione, &one, c_pntr,
                 &ione);
      }

      val_pntr += iblk_size;
    }
  }

} /* dvbr_sparax_basic */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void AZ_VBR_matvec_mult(double *b, double *c, AZ_MATRIX *Amat,int proc_config[])



/******************************************************************************

  c = Ab:
  Sparse (square) overlapped matrix-vector multiply, using the distributed
  variable block row (VBR) data structure.

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============


  b:               Contains the vector b.

  c:               Contains the result vector c.



  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  Amat:            Structure used to represent the matrix (see file az_aztec.h
                   and Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

******************************************************************************/

{

  /* local variables */

  int          num_blks;
  int exchange_flag;
  double *val;
  int *data_org, *bindx, *rpntr,*cpntr,*bpntr;
  val = Amat->val;
  bindx = Amat->bindx;
  rpntr =Amat->rpntr;
  cpntr = Amat->cpntr;
  bpntr = Amat->bpntr;
  data_org = Amat->data_org;
  exchange_flag = 1;



  /**************************** execution begins *****************************/

  num_blks = data_org[AZ_N_int_blk];

  /* multiple processors */

  if (exchange_flag )
     num_blks += data_org[AZ_N_bord_blk];

  /*
   * It is possible to overlap the communication and computation here for
   * performace gains.  The idea is to gather the messages each processor
   * needs to send (AZ_gather_mesg_info), write these out (send them),
   * perform a 'dvbr_sparax_basic' on the INTERNAL portion of the sparse
   * matrix-vector product, read the messages from the neighboring processors
   * and then perform a 'dvbr_sparax_basic' on the BOUNDARY portion of the
   * product.  We do not support this capability at this time.  SAH, 2/1996
   */

  /* perform the sparax - NOTE: the boundary exchange is done inside
   dvbr_sparax_basic */

  dvbr_sparax_basic(num_blks, val, bindx, rpntr, cpntr, bpntr, b, c,
        exchange_flag, data_org, proc_config);

}
/* AZ_VBR_matvec_mult */


void AZ_MSR_matvec_mult (double *b, double *c,AZ_MATRIX *Amat,int proc_config[])


/******************************************************************************
  c = Ab:
  Sparse (square) overlapped matrix-vector multiply, using the  MSR
  data structure .

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============


  b:               Contains the vector b.

  c:               Contains the result vector c.



  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  Amat:            Structure used to represent the matrix (see file az_aztec.h
                   and Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

******************************************************************************/



{
  double *val;
  int *data_org, *bindx;
  register int j, irow;
 int          N;
#ifndef AZ_DONT_UNROLL_LOOPS
 int          *bindx_ptr;
 double       sum;
#else
 int nzeros, bindx_row, k;
#endif

  val = Amat->val;
  bindx = Amat->bindx;
  data_org = Amat->data_org;


  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  /* exchange boundary info */

  AZ_exchange_bdry(b, data_org, proc_config);

  /* This is the default */
#ifndef AZ_DONT_UNROLL_LOOPS

  j = bindx[0];
  bindx_ptr = &bindx[j];
  for (irow = 0; irow < N; irow++) {
    sum =  val[irow]*b[irow];
    while (j+10 < bindx[irow+1]) {
      sum += val[j+9]*b[bindx_ptr[9]] +
	val[j+8]*b[bindx_ptr[8]] +
	val[j+7]*b[bindx_ptr[7]] +
	val[j+6]*b[bindx_ptr[6]] +
	val[j+5]*b[bindx_ptr[5]] +
	val[j+4]*b[bindx_ptr[4]] +
	val[j+3]*b[bindx_ptr[3]] +
	val[j+2]*b[bindx_ptr[2]] +
	val[j+1]*b[bindx_ptr[1]] +
	val[j]*b[*bindx_ptr];
      bindx_ptr += 10;
      j += 10;
    }
    while (j < bindx[irow+1]) {
      sum += val[j++] * b[*bindx_ptr++];
    }
    c[irow] = sum;
  }

  /* This is available for backward compatibility.  Turn on by specifying -DAZ_DONT_UNROLL_LOOPS to the compiler */
#else                                                                           

  for (irow = 0; irow < N; irow++) {
    /* compute diagonal contribution */
    *c = val[irow] * b[irow];
    /* nonzero off diagonal contribution */
    bindx_row = bindx[irow];
    nzeros    = bindx[irow+1] - bindx_row;
    for (j = 0; j < nzeros; j++) {
      k   = bindx_row + j;
      *c += val[k] * b[bindx[k]];
    }
    c++;
    
  }
  
#endif

} /* AZ_MSR_matvec_mult */

