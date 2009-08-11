/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
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

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_SingletonFilter.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"

//==============================================================================
Tifpack_SingletonFilter::Tifpack_SingletonFilter(const Teuchos::RefCountPtr<Tpetra_RowMatrix>& Matrix) :
  A_(Matrix),
  NumSingletons_(0),
  NumRows_(0),
  NumRowsA_(0),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0),
  NumNonzeros_(0)
{
  // use this filter only on serial matrices
  if (A_->Comm().NumProc() != 1) {
    cerr << "Tifpack_SingletonFilter can be used with Comm().NumProc() == 1" << endl;
    cerr << "only. This class is a tool for Tifpack_AdditiveSchwarz," << endl;
    cerr << "and it is not meant to be used otherwise." << endl;
    exit(EXIT_FAILURE);
  }
  
  if ((A_->NumMyRows() != A_->NumGlobalRows()) ||
     (A_->NumMyRows() != A_->NumMyCols()))
    TIFPACK_CHK_ERRV(-1);
  
  NumRowsA_ = (A_->NumMyRows());
  MaxNumEntriesA_ = A_->MaxNumEntries();

  Indices_.resize(MaxNumEntriesA_);
  Values_.resize(MaxNumEntriesA_);
  Reorder_.resize(A_->NumMyRows());

  for (int i = 0 ; i < NumRowsA_ ; ++i)
    Reorder_[i] = -1;

  // first check how may singletons I do have
  for (int i = 0 ; i < NumRowsA_ ; ++i) {
    int Nnz;
    TIFPACK_CHK_ERRV(A_->ExtractMyRowCopy(i,MaxNumEntriesA_,Nnz,
					&Values_[0], &Indices_[0]));
    if (Nnz != 1) {
      Reorder_[i] = NumRows_++;
    }
    else {
      NumSingletons_++;
    }
  }

  InvReorder_.resize(NumRows_);
  for (int i = 0 ; i < NumRowsA_ ; ++i) {
    if (Reorder_[i] < 0)
      continue;
    InvReorder_[Reorder_[i]] = i;
  }
 
  NumEntries_.resize(NumRows_);
  SingletonIndex_.resize(NumSingletons_);

  // now compute the nonzeros per row
  int count = 0;
  for (int i = 0 ; i < A_->NumMyRows() ; ++i) {

    int Nnz;
    TIFPACK_CHK_ERRV(A_->ExtractMyRowCopy(i,MaxNumEntriesA_,Nnz,
					  &Values_[0], &Indices_[0]));
    int ii = Reorder_[i];
    if (ii >= 0) {
      assert (Nnz != 1);

      NumEntries_[ii] = Nnz;
      NumNonzeros_ += Nnz;
      if (Nnz > MaxNumEntries_)
	MaxNumEntries_ = Nnz;
    }
    else {
      SingletonIndex_[count] = i;
      count++;
    }
  }

  Map_ = Teuchos::rcp( new Tpetra_Map(NumRows_,0,Comm()) );

  // and finish up with the diagonal entry
  Diagonal_ = Teuchos::rcp( new Tpetra_Vector(*Map_) );

  Tpetra_Vector Diagonal(A_->Map());
  A_->ExtractDiagonalCopy(Diagonal);
  for (int i = 0 ; i < NumRows_ ; ++i) {
    int ii = InvReorder_[i];
    (*Diagonal_)[i] = Diagonal[ii];
  }
  
}

//==============================================================================
int Tifpack_SingletonFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{
  int Nnz;

  if (Length < NumEntries_[MyRow])
    TIFPACK_CHK_ERR(-1);

  int Row = InvReorder_[MyRow];
  TIFPACK_CHK_ERR(A_->ExtractMyRowCopy(Row,MaxNumEntriesA_,Nnz,
				     &Values_[0],&Indices_[0]));
  NumEntries = 0;
  for (int i = 0 ; i < Nnz ; ++i) {
    int ii = Reorder_[Indices_[i]];
    if ( ii >= 0) {
      Indices[NumEntries] = ii;
      Values[NumEntries] = Values_[i];
      NumEntries++;
    }
  }
  return(0);   
}

//==============================================================================
int Tifpack_SingletonFilter::
ExtractDiagonalCopy(Tpetra_Vector & Diagonal) const
{
  Diagonal = *Diagonal_;
  return(0);
}

//==============================================================================
int Tifpack_SingletonFilter::
Multiply(bool TransA, const Tpetra_MultiVector& X, 
	 Tpetra_MultiVector& Y) const
{

  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors())
    TIFPACK_CHK_ERR(-1);

  Y.PutScalar(0.0);

  vector<int> Indices(MaxNumEntries_);
  vector<double> Values(MaxNumEntries_);

  // cycle over all rows of the original matrix
  for (int i = 0 ; i < A_->NumMyRows() ; ++i) {

    if (Reorder_[i] < 0)
      continue; // skip singleton rows
    
    int Nnz;
    A_->ExtractMyRowCopy(i,MaxNumEntriesA_,Nnz,
		     &Values[0], &Indices[0]);
    if (!TransA) {
      // no transpose first
      for (int j = 0 ; j < NumVectors ; ++j) {
	for (int k = 0 ; k < Nnz ; ++k) {
	  if (Reorder_[i] >= 0)
	    Y[j][i] += Values[k] * X[j][Reorder_[Indices[k]]];
	}
      }
    }
    else {
      // transpose here
      for (int j = 0 ; j < NumVectors ; ++j) {
	for (int k = 0 ; k < Nnz ; ++k) {
	  if (Reorder_[i] >= 0)
	    Y[j][Reorder_[Indices[k]]] += Values[k] * X[j][i];
	}
      }
    }
  }

  return(0);
}

//==============================================================================
int Tifpack_SingletonFilter::
Solve(bool Upper, bool Trans, bool UnitDiagonal, 
      const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  return(-1);
}

//==============================================================================
int Tifpack_SingletonFilter::
Apply(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  TIFPACK_CHK_ERR(Multiply(false,X,Y));
  return(0);
}

//==============================================================================
int Tifpack_SingletonFilter::
ApplyInverse(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const
{
  return(-1); // NOT IMPLEMENTED AT THIS STAGE
}

//==============================================================================
int Tifpack_SingletonFilter::
SolveSingletons(const Tpetra_MultiVector& RHS, 
		Tpetra_MultiVector& LHS)
{
  for (int i = 0 ; i < NumSingletons_ ; ++i) {
    int ii = SingletonIndex_[i];
    // get the diagonal value for the singleton
    int Nnz;
    A_->ExtractMyRowCopy(ii,MaxNumEntriesA_,Nnz,
			&Values_[0], &Indices_[0]);
    for (int j = 0 ; j < Nnz ; ++j) {
      if (Indices_[j] == ii) {
	for (int k = 0 ; k < LHS.NumVectors() ; ++k)
	  LHS[k][ii] = RHS[k][ii] / Values_[j];
      }
    }
  }

  return(0);
}

//==============================================================================
int Tifpack_SingletonFilter::
CreateReducedRHS(const Tpetra_MultiVector& LHS,
		 const Tpetra_MultiVector& RHS, 
		 Tpetra_MultiVector& ReducedRHS)
{
  int NumVectors = LHS.NumVectors();

  for (int i = 0 ; i < NumRows_ ; ++i)
    for (int k = 0 ; k < NumVectors ; ++k)
      ReducedRHS[k][i] = RHS[k][InvReorder_[i]];

  for (int i = 0 ; i < NumRows_ ; ++i) {
    int ii = InvReorder_[i];
    int Nnz;
    TIFPACK_CHK_ERR(A_->ExtractMyRowCopy(ii,MaxNumEntriesA_,Nnz,
					&Values_[0], &Indices_[0]));

    for (int j = 0 ; j < Nnz ; ++j) {
      if (Reorder_[Indices_[j]] == -1) {
	for (int k = 0 ; k < NumVectors ; ++k)
	  ReducedRHS[k][i] -= Values_[j] * LHS[k][Indices_[j]];
      }
    }
  }
  return(0);
}

//==============================================================================
int Tifpack_SingletonFilter::
UpdateLHS(const Tpetra_MultiVector& ReducedLHS,
	  Tpetra_MultiVector& LHS)
{
  for (int i = 0 ; i < NumRows_ ; ++i)
    for (int k = 0 ; k < LHS.NumVectors() ; ++k)
      LHS[k][InvReorder_[i]] = ReducedLHS[k][i];

  return(0);
}
