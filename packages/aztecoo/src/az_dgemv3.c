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
