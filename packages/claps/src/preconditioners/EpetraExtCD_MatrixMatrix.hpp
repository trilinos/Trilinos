//@HEADER
// ************************************************************************
//
//         Claps: A Collection of Domain Decomposition Preconditioners
//                and Solvers
//         Copyright (2006) Sandia Corporation
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
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef EPETRAEXTCD_MATRIXMATRIX_H
#define EPETRAEXTCD_MATRIXMATRIX_H

#include <Epetra_RowMatrix.h>
#include <Epetra_CrsMatrix.h>

namespace EpetraExtCD {

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
    @param B Input, must already have had 'FillComplete()' called.
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
			const Epetra_CrsMatrix& B,
			Epetra_CrsMatrix* & C);

  };//class MatrixMatrix
}//namespace EpetraExtCD

#endif

