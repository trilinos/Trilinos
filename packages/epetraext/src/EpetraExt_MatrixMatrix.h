// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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
// @HEADER
#ifndef EPETRAEXT_MATRIXMATRIX_H
#define EPETRAEXT_MATRIXMATRIX_H

#include <EpetraExt_ConfigDefs.h>

class Epetra_CrsMatrix;

namespace EpetraExt {

  /** Collection of matrix-matrix operations. This class basically
      functions as a namespace, containing only static methods.
      See the program epetraext/test/MatrixMatrix/cxx_main.cpp for
      a usage example.
   */
class MatrixMatrix {

  public:
    /** destructor */
    virtual ~MatrixMatrix(){}

    /** Given Epetra_CrsMatrix objects A, B and C, form the product C = A*B.
	In a parallel setting, A and B need not have matching distributions,
	but C needs to have the same row-map as A.

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. On entry to this method, it doesn't matter whether
             FillComplete() has already been called on C or not. If it has,
	     then C's graph must already contain all nonzero locations that
	     will be produced when forming the product A*B. On exit,
	     C.FillComplete() will have been called.

    @return error-code, 0 if successful. non-zero returns may result if A or
             B are not already Filled, or if errors occur in putting values
             into C, etc.
     */
    static int Multiply(const Epetra_CrsMatrix& A,
			bool transposeA,
			const Epetra_CrsMatrix& B,
			bool transposeB,
			Epetra_CrsMatrix& C);

    /** Given Epetra_CrsMatrix objects A and B, form the sum B = a*A + b*B

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param scalarA Input, scalar multiplier for matrix A.
    @param B Result. On entry to this method, it doesn't matter whether
             FillComplete() has already been called on B or not. If it has,
	     then B's graph must already contain all nonzero locations that
	     will be produced when forming the sum.
    @param scalarB Input, scalar multiplier for matrix B.

    @return error-code, 0 if successful. non-zero returns may result if A is
             not already Filled, or if errors occur in putting values
             into B, etc.
     */
    static int Add(const Epetra_CrsMatrix& A,
                   bool transposeA,
                   double scalarA,
                   Epetra_CrsMatrix& B,
                   double scalarB);

};//class MatrixMatrix


/**
 *Method for internal use... sparsedot forms a dot-product between two
 *sparsely-populated 'vectors'.
 *Important assumption: assumes the indices in u_ind and v_ind are sorted.
 */
 double sparsedot(double* u, int* u_ind, int u_len,
		  double* v, int* v_ind, int v_len);
}//namespace EpetraExt

#endif

