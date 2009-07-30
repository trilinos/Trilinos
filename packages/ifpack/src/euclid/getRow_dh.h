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

#ifndef GET_ROW_DH
#define GET_ROW_DH

#include "euclid_common.h"
#include "call_epetra.h"

/* "row" refers to global row number */

#ifdef __cplusplus
extern "C" {
#endif

extern void EuclidGetDimensions(void *A, int *beg_row, int *rowsLocal, int *rowsGlobal);
extern void EuclidGetRow(void *A, int row, int *len, int **ind, double **val);
extern void EuclidRestoreRow(void *A, int row, int *len, int **ind, double **val);

extern int EuclidReadLocalNz(void *A);

extern void PrintMatUsingGetRow(void* A, int beg_row, int m,
                          int *n2o_row, int *n2o_col, char *filename);

#ifdef __cplusplus
}
#endif
#endif

