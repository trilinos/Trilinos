/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_ReorderFilter.h"
#include "Ifpack_Reordering.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//==============================================================================
Ifpack_ReorderFilter::Ifpack_ReorderFilter(const Teuchos::RefCountPtr<Epetra_RowMatrix>& Matrix_in,
                                           const Teuchos::RefCountPtr<Ifpack_Reordering>& Reordering_in) :
  A_(Matrix_in),
  Reordering_(Reordering_in),
  NumMyRows_(Matrix_in->NumMyRows()),
  MaxNumEntries_(Matrix_in->MaxNumEntries())
{
}

//==============================================================================
Ifpack_ReorderFilter::Ifpack_ReorderFilter(const Ifpack_ReorderFilter& RHS) :
  A_(Matrix()),
  Reordering_(Reordering()),
  NumMyRows_(RHS.NumMyRows()),
  MaxNumEntries_(RHS.MaxNumEntries())
{
  strcpy(Label_,RHS.Label());
}

//==============================================================================
Ifpack_ReorderFilter& 
Ifpack_ReorderFilter::operator=(const Ifpack_ReorderFilter& RHS)
{
  if (this == &RHS)
    return (*this);

  A_ = RHS.Matrix();

  Reordering_ = RHS.Reordering();
  MaxNumEntries_ = RHS.MaxNumEntries();
  NumMyRows_ = RHS.NumMyRows();

  strcpy(Label_,RHS.Label());
  return(*this);
}

//==============================================================================
int Ifpack_ReorderFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
                 double *Values, int * Indices) const
{
  int MyReorderdRow = Reordering_->InvReorder(MyRow);

  IFPACK_CHK_ERR(Matrix()->ExtractMyRowCopy(MyReorderdRow,MaxNumEntries_,
                                           NumEntries, Values,Indices));

  // suppose all elements are local. Note that now
  // Indices can have indices in non-increasing order.
  for (int i = 0 ; i < NumEntries ; ++i) {
    Indices[i] = Reordering_->Reorder(Indices[i]);
  }

  return(0);
}

//==============================================================================
int Ifpack_ReorderFilter::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  Epetra_Vector DiagonalTilde(Diagonal.Map());
  IFPACK_CHK_ERR(Matrix()->ExtractDiagonalCopy(DiagonalTilde));
  IFPACK_CHK_ERR((Reordering_->P(DiagonalTilde,Diagonal)));
  return(0);
}

//==============================================================================
int Ifpack_ReorderFilter::
Multiply(bool TransA, const Epetra_MultiVector& X, 
         Epetra_MultiVector& Y) const
{
  // need two additional vectors
  Epetra_MultiVector Xtilde(X.Map(),X.NumVectors());
  Epetra_MultiVector Ytilde(Y.Map(),Y.NumVectors());
  // bring X back to original ordering
  Reordering_->Pinv(X,Xtilde);
  // apply original matrix
  IFPACK_CHK_ERR(Matrix()->Multiply(TransA,Xtilde,Ytilde));
  // now reorder result
  Reordering_->P(Ytilde,Y);


  return(0);
}

//==============================================================================
int Ifpack_ReorderFilter::
Solve(bool Upper, bool Trans, bool UnitDiagonal, 
      const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-98);
}

//==============================================================================
int Ifpack_ReorderFilter::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int ierr = Multiply(UseTranspose(),X,Y);
  IFPACK_RETURN(ierr);
}
