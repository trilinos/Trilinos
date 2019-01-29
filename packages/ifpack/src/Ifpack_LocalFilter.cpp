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

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Ifpack_LocalFilter.h"

//==============================================================================
Ifpack_LocalFilter::Ifpack_LocalFilter(const Teuchos::RefCountPtr<const Epetra_RowMatrix>& Matrix) :
  Matrix_(Matrix),
  NumRows_(0),
  NumNonzeros_(0),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0)
{
  sprintf(Label_,"%s","Ifpack_LocalFilter");

#ifdef HAVE_MPI
  SerialComm_ = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_SELF) );
#else
  SerialComm_ = Teuchos::rcp( new Epetra_SerialComm );
#endif

  // localized matrix has all the local rows of Matrix
  NumRows_ = Matrix->NumMyRows();

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  // build a linear map, based on the serial communicator
  Map_ = Teuchos::rcp( new Epetra_Map(NumRows_,0,*SerialComm_) );
#endif

  // NumEntries_ will contain the actual number of nonzeros
  // for each localized row (that is, without external nodes,
  // and always with the diagonal entry)
  NumEntries_.resize(NumRows_);

  // want to store the diagonal vector. FIXME: am I really useful?
  Diagonal_ = Teuchos::rcp( new Epetra_Vector(*Map_) );
  if (Diagonal_ == Teuchos::null) IFPACK_CHK_ERRV(-5);
  
  // store this for future access to ExtractMyRowCopy().
  // This is the # of nonzeros in the non-local matrix
  MaxNumEntriesA_ = Matrix->MaxNumEntries();
  // tentative value for MaxNumEntries. This is the number of
  // nonzeros in the local matrix
  MaxNumEntries_ = Matrix->MaxNumEntries();

  // ExtractMyRowCopy() will use these vectors
  Indices_.resize(MaxNumEntries_);
  Values_.resize(MaxNumEntries_);

  // now compute:
  // - the number of nonzero per row
  // - the total number of nonzeros
  // - the diagonal entries

  // compute nonzeros (total and per-row), and store the
  // diagonal entries (already modified)
  int ActualMaxNumEntries = 0;

  for (int i = 0 ; i < NumRows_ ; ++i) {
    
    NumEntries_[i] = 0;
    int Nnz, NewNnz = 0;
    IFPACK_CHK_ERRV(ExtractMyRowCopy(i,MaxNumEntries_,Nnz,&Values_[0],&Indices_[0]));

    for (int j = 0 ; j < Nnz ; ++j) {
      if (Indices_[j] < NumRows_ ) ++NewNnz;

      if (Indices_[j] == i)
	(*Diagonal_)[i] = Values_[j];
    }

    if (NewNnz > ActualMaxNumEntries)
      ActualMaxNumEntries = NewNnz;

    NumNonzeros_ += NewNnz;
    NumEntries_[i] = NewNnz;

  }
 
  MaxNumEntries_ = ActualMaxNumEntries;
}

//==============================================================================
int Ifpack_LocalFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{
  if ((MyRow < 0) || (MyRow >= NumRows_)) {
    IFPACK_CHK_ERR(-1); // range not valid
  }

  if (Length < NumEntries_[MyRow])
    return(-1);

  // always extract using the object Values_ and Indices_.
  // This is because I need more space than that given by
  // the user (for the external nodes)
  int Nnz;
  int ierr = Matrix_->ExtractMyRowCopy(MyRow,MaxNumEntriesA_,Nnz,
				       &Values_[0],&Indices_[0]);

  IFPACK_CHK_ERR(ierr);

  // populate the user's vectors, add diagonal if not found
  NumEntries = 0;

  for (int j = 0 ; j < Nnz ; ++j) {
    // only local indices
    if (Indices_[j] < NumRows_ ) {
      Indices[NumEntries] = Indices_[j];
      Values[NumEntries] = Values_[j];
      ++NumEntries;
    }
  }
    
  return(0);

}

//==============================================================================
int Ifpack_LocalFilter::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  if (!Diagonal.Map().SameAs(*Map_))
    IFPACK_CHK_ERR(-1);
  Diagonal = *Diagonal_;
  return(0);
}

//==============================================================================
int Ifpack_LocalFilter::Apply(const Epetra_MultiVector& X,
	  Epetra_MultiVector& Y) const 
{

  // skip expensive checks, I suppose input data are ok

  Y.PutScalar(0.0);
  int NumVectors = Y.NumVectors();

  double** X_ptr;
  double** Y_ptr;
  X.ExtractView(&X_ptr);
  Y.ExtractView(&Y_ptr);

  for (int i = 0 ; i < NumRows_ ; ++i) {
    
    int Nnz;
    int ierr = Matrix_->ExtractMyRowCopy(i,MaxNumEntriesA_,Nnz,&Values_[0],
                                         &Indices_[0]);
    IFPACK_CHK_ERR(ierr);

    for (int j = 0 ; j < Nnz ; ++j) {
      if (Indices_[j] < NumRows_ ) {
	for (int k = 0 ; k < NumVectors ; ++k)
	  Y_ptr[k][i] += Values_[j] * X_ptr[k][Indices_[j]];
      }
    }
  }

  return(0);
}

//==============================================================================
int Ifpack_LocalFilter::ApplyInverse(const Epetra_MultiVector& /* X */,
		 Epetra_MultiVector& /* Y */) const
{
  IFPACK_CHK_ERR(-1); // not implemented
}

//==============================================================================
const Epetra_BlockMap& Ifpack_LocalFilter::Map() const
{
  return(*Map_);
}
