//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
#ifndef EPETRAEXT_MATRIXMATRIX_H
#define EPETRAEXT_MATRIXMATRIX_H

#include <EpetraExt_ConfigDefs.h>

class Epetra_CrsMatrix;
class Epetra_Map;

#ifdef HAVE_VECTOR
#include <vector>
#endif

namespace EpetraExt {
  class CrsMatrixStruct;


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
	     C.FillComplete() will have been called, unless the last argument
             to this function is specified to be false.
    @param call_FillComplete_on_result Optional argument, defaults to true.
           Power users may specify this argument to be false if they *DON'T*
           want this function to call C.FillComplete. (It is often useful
           to allow this function to call C.FillComplete, in cases where
           one or both of the input matrices are rectangular and it is not
           trivial to know which maps to use for the domain- and range-maps.)

    @return error-code, 0 if successful. non-zero returns may result if A or
             B are not already Filled, or if errors occur in putting values
             into C, etc.
     */
    static int Multiply(const Epetra_CrsMatrix& A,
			bool transposeA,
			const Epetra_CrsMatrix& B,
			bool transposeB,
			Epetra_CrsMatrix& C,
                        bool call_FillComplete_on_result=true);

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

    /** Given Epetra_CrsMatrix objects A and B, form the sum C = a*A + b*B

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param scalarA Input, scalar multiplier for matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param scalarB Input, scalar multiplier for matrix B.
    @param C Result. On entry to this method, C can be NULL or a pointer
             to an unfilled or filled matrix. If C is NULL then a new
             object is allocated and must be deleted by the user.
             If C is not NULL and FillComplete has already
             been called then the sparsity pattern is assumed to be fixed
             and compatible  with the sparsity of A+B. If FillComplete has
             not been called then the sum is completed and the function
             returns without calling FillComplete on C.

    @return error-code, 0 if successful. non-zero returns may result if A or is
             not already Filled, or if errors occur in putting values
             into C, etc.
     */
    static int Add(const Epetra_CrsMatrix& A,
                   bool transposeA,
                   double scalarA,
                   const Epetra_CrsMatrix & B,
                   bool transposeB,
                   double scalarB,
                   Epetra_CrsMatrix * & C);

 private:
    static int mult_A_B_reuse(const Epetra_CrsMatrix & A,
		       const Epetra_CrsMatrix & B,
		       CrsMatrixStruct& Bview,
		       std::vector<int> & Bcol2Ccol,
		       std::vector<int> & Bimportcol2Ccol,
		       Epetra_CrsMatrix& C);
    
    static int mult_A_B_newmatrix(const Epetra_CrsMatrix & A,
			   const Epetra_CrsMatrix & B,
			   CrsMatrixStruct& Bview,
			   std::vector<int> & Bcol2Ccol,
			   std::vector<int> & Bimportcol2Ccol,
			   Epetra_CrsMatrix& C);

    static int mult_A_B(const Epetra_CrsMatrix & A,
		 CrsMatrixStruct & Aview,
		 const Epetra_CrsMatrix & B,
		 CrsMatrixStruct& Bview,
		 Epetra_CrsMatrix& C,
		 bool call_FillComplete_on_result);

    static int ReplaceCrsGraphData(Epetra_CrsMatrix & C, const Epetra_Map & DomainMap, const Epetra_Map & RangeMap);

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

