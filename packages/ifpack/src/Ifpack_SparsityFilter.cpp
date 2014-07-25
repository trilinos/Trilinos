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
#include "Ifpack_SparsityFilter.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//==============================================================================
Ifpack_SparsityFilter::Ifpack_SparsityFilter(const Teuchos::RefCountPtr<Epetra_RowMatrix>& Matrix,
					     int AllowedEntries, 
					     int AllowedBandwidth) :
  A_(Matrix),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0),
  AllowedBandwidth_(AllowedBandwidth),
  AllowedEntries_(AllowedEntries),
  NumNonzeros_(0),
  NumRows_(0)
{
  // use this filter only on serial matrices
  if (A_->Comm().NumProc() != 1) {
    cerr << "Ifpack_SparsityFilter can be used with Comm().NumProc() == 1" << endl;
    cerr << "only. This class is a tool for Ifpack_AdditiveSchwarz," << endl;
    cerr << "and it is not meant to be used otherwise." << endl;
    exit(EXIT_FAILURE);
  }

  // only square serial matrices
  if ((A_->NumMyRows() != A_->NumMyCols()) ||
     (A_->NumMyRows() != A_->NumGlobalRows64()))
    IFPACK_CHK_ERRV(-1);

  NumRows_ = A_->NumMyRows();
  MaxNumEntriesA_ = A_->MaxNumEntries();
  Indices_.resize(MaxNumEntriesA_);
  Values_.resize(MaxNumEntriesA_);

  // default value is to not consider bandwidth
  if (AllowedBandwidth_ == -1)
    AllowedBandwidth_ = NumRows_;
  
  // computes the number of nonzero elements per row in the 
  // dropped matrix. Stores this number in NumEntries_.
  // Also, computes the global number of nonzeros.
  std::vector<int>    Ind(MaxNumEntriesA_);
  std::vector<double> Val(MaxNumEntriesA_);

  NumEntries_.resize(NumRows_);
  for (int i = 0 ; i < NumRows_ ; ++i)
    NumEntries_[i] = MaxNumEntriesA_;

  for (int i = 0 ; i < A_->NumMyRows() ; ++i) {
    int Nnz;
    IFPACK_CHK_ERRV(ExtractMyRowCopy(i,MaxNumEntriesA_,Nnz,
				     &Val[0], &Ind[0]));

    NumEntries_[i] = Nnz;
    NumNonzeros_ += Nnz;
    if (Nnz > MaxNumEntries_)
      MaxNumEntries_ = Nnz;

  }
}

//==============================================================================
int Ifpack_SparsityFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{
  if (Length < NumEntries_[MyRow])
    IFPACK_CHK_ERR(-1);

  int Nnz;
  IFPACK_CHK_ERR(A_->ExtractMyRowCopy(MyRow,MaxNumEntriesA_,Nnz,
				     &Values_[0],&Indices_[0]));

  double Threshold = 0.0;
    
  // this `if' is to define the cut-off value
  if (Nnz > AllowedEntries_) {
 
    std::vector<double> Values2(Nnz);
    int count = 0;
    for (int i = 0 ; i < Nnz ; ++i) {
      // skip diagonal entry (which is always inserted)
      if (Indices_[i] == MyRow)
	continue;
      // put absolute value
      Values2[count] = IFPACK_ABS(Values_[i]);
      count++;
    }

    for (int i = count ; i < Nnz ; ++i)
      Values2[i] = 0.0;

    // sort in descending order
    std::sort(Values2.rbegin(),Values2.rend());
    // get the cut-off value
    Threshold = Values2[AllowedEntries_ - 1];

  }

  // loop over all nonzero elements of row MyRow,
  // and drop elements below specified threshold.
  // Also, add value to zero diagonal
  NumEntries = 0;

  for (int i = 0 ; i < Nnz ; ++i) {

    if (IFPACK_ABS(Indices_[i] - MyRow) > AllowedBandwidth_)
      continue;

    if ((Indices_[i] != MyRow) && (IFPACK_ABS(Values_[i]) < Threshold))
      continue;

    Values[NumEntries] = Values_[i];
    Indices[NumEntries] = Indices_[i];

    NumEntries++;
    if (NumEntries > AllowedEntries_)
      break;
  }

  return(0);
}

//==============================================================================
int Ifpack_SparsityFilter::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  int ierr = A_->ExtractDiagonalCopy(Diagonal);
  IFPACK_RETURN(ierr);
}

//==============================================================================
int Ifpack_SparsityFilter::
Multiply(bool TransA, const Epetra_MultiVector& X, 
	 Epetra_MultiVector& Y) const
{

  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors())
    IFPACK_CHK_ERR(-1);

  Y.PutScalar(0.0);

  std::vector<int> Indices(MaxNumEntries_);
  std::vector<double> Values(MaxNumEntries_);

  for (int i = 0 ; i < A_->NumMyRows() ; ++i) {

    int Nnz;
    ExtractMyRowCopy(i,MaxNumEntries_,Nnz,
		     &Values[0], &Indices[0]);
    if (!TransA) {
      // no transpose first
      for (int j = 0 ; j < NumVectors ; ++j) {
	for (int k = 0 ; k < Nnz ; ++k) {
	  Y[j][i] += Values[k] * X[j][Indices[k]];
	}
      }
    }
    else {
      // transpose here
      for (int j = 0 ; j < NumVectors ; ++j) {
	for (int k = 0 ; k < Nnz ; ++k) {
	  Y[j][Indices[k]] += Values[k] * X[j][i];
	}
      }
    }
  }

  return(0);
}

//==============================================================================
int Ifpack_SparsityFilter::
Solve(bool Upper, bool Trans, bool UnitDiagonal, 
      const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-98);
}

//==============================================================================
int Ifpack_SparsityFilter::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int ierr = Multiply(UseTranspose(),X,Y);
  IFPACK_RETURN(ierr);
}

//==============================================================================
int Ifpack_SparsityFilter::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(-98); 
}
