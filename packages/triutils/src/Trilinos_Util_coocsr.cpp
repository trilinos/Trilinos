//@HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2001) Sandia Corporation
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

#include "Trilinos_Util.h"

int Trilinos_Util_coocsr( int nrow, int nnz, double *a, int *ir, int *jc,
                          double *ao, int *jao, int *iao) {
    int i, j, k;
    double x;
    int k0, iad;
/*
 -----------------------------------------------------------------------
  Coordinate     to   Compressed Sparse Row 
 -----------------------------------------------------------------------
 
 converts a matrix that is stored in coordinate format 
  a, ir, jc into a row general sparse ao, jao, iao format. 

 on entry: 
 --------- 
 nrow	= dimension of the matrix 
 nnz	= number of nonzero elements in matrix 
 a, 
 ir, 
 jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz 

         nonzero elements of the matrix with a(k) = actual real value of
 
 	  the elements, ir(k) = its row number and jc(k) = its column 
 	  number. The order of the elements is arbitrary. 

 on return: 
 ----------- 
 ir 	is destroyed 

 ao, jao, iao = matrix in general sparse matrix format with ao 
 	continung the real values, jao containing the column indices, 
 	and iao being the pointer to the beginning of the row, 
 	in arrays ao, jao. 

 Notes: 
 ------ This routine is NOT in place.  See coicsr 

----------------------------------------------------------------------- -
*/
    for (k = 0; k <= nrow; k++) iao[k] = 0;

 // determine row-lengths. 

    for (k = 0; k <nnz ; k++) ++iao[ir[k]];

 // starting position of each row.. 

    k = 0;
    for (j = 0; j <= nrow; j++) {
	k0 = iao[j];
	iao[j] = k;
	k += k0;
    }
 
 // go through the structure  once more. Fill in output matrix. 

    for (k = 0; k <nnz; k++) {
	i = ir[k];
	j = jc[k];
	x = a[k];
	iad = iao[i];
	ao[iad] = x;
	jao[iad] = j;
	iao[i] = iad + 1;
    }
 
  // shift back iao 
    for (j = nrow-1; j >= 0; j--) {
	iao[j + 1] = iao[j];
    }
    iao[0] = 0;
    return 0;
}

