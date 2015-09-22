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

extern void AZ_dgemv3(int, int, double *, double *, double *);


void AZ_dgemv3(int m, int n, double *a, double *x, double *y)

/*******************************************************************************

  Perform the matrix-vector operation y = y - Ax where x and y are vectors and
  A is an m-by-n matrix.  This modeled on the blas 2 routine dgemv and used for
  small matrix sizes.  Also, there is a special case for m = n = 5, (3d fluids).

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  m:               Number of rows of the matrix a.

  n:               Number of columns of the matrix a.

  a:               The leading m-by-n part of the array a must contain the
                   matrix of coefficients.

  x:               The vector of length at least n to multiply by A.

  y:               The vector of length at least m which, on output, will
                   contain the result of the -Ax operation summed into it.
                   NOTE:  the result is the SUM of the entry y and -Ax.

*******************************************************************************/

{

  /* local variables */

  register int i = 0;

  /**************************** execution begins ******************************/

  switch (m) {
  case 5:
    while (i++ < n) {
      *y     -= *a++ * *x;
      *(y+1) -= *a++ * *x;
      *(y+2) -= *a++ * *x;
      *(y+3) -= *a++ * *x;
      *(y+4) -= *a++ * *x;

      x++;

    }

    break;

  default:

     while (i++ < n) {
      register int j = 0;
      while (j++ < m)
        *y++ -= *a++ * *x;
      y -= m; x++;
    }
  }

} /* AZ_dgemv3 */
