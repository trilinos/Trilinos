/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

/* shell sort adopted from Edmond Chow */

#include "shellSort_dh.h"

#undef __FUNC__
#define __FUNC__ "shellSort_int"
void shellSort_int(const int n, int *x)
{
  START_FUNC_DH
  int m, max, j, k, itemp;

  m = n/2;
  while (m > 0) {
    max = n - m;
    for (j=0; j<max; j++) {
      for (k=j; k>=0; k-=m) {
        if (x[k+m] >= x[k]) break;
        itemp = x[k+m];
        x[k+m] = x[k];
        x[k] = itemp;
      }
    }
    m = m/2;
  }
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "shellSort_float"
void shellSort_float(const int n, double *x)
{
  START_FUNC_DH
  int m, max, j, k;
  double itemp;

  m = n/2;
  while (m > 0) {
    max = n - m;
    for (j=0; j<max; j++) {
      for (k=j; k>=0; k-=m) {
        if (x[k+m] >= x[k]) break;
        itemp = x[k+m];
        x[k+m] = x[k];
        x[k] = itemp;
      }
    }
    m = m/2;
  }
  END_FUNC_DH
}


#if 0
#undef __FUNC__
#define __FUNC__ "shellSort_int_float"
void shellSort_int_float(int n, int *x, VAL_DH *xVals)
{
  START_FUNC_DH
  int m, max, j, k, itemp;
  VAL_DH atemp;

  m = n/2;
  while (m > 0) {
    max = n - m;
    for (j=0; j<max; j++) {
      for (k=j; k>=0; k-=m) {
        if (x[k+m] >= x[k]) break;
        itemp = x[k+m];
        atemp = xVals[k+m];
        x[k+m] = x[k];
        /* xVals[k+m] = xVals[k]; */
        x[k] = itemp;
        xVals[k] = atemp;
      }
    }
    m = m/2;
  }
  END_FUNC_DH
}
#endif
