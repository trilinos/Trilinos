/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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
  IFPACK_RETURN(Multiply(UseTranspose(),X,Y));
}
