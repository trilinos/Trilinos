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

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_BlockMap.hpp"
#include "Tifpack_LocalFilter.hpp"

//==============================================================================
Tifpack_LocalFilter::Tifpack_LocalFilter(const Teuchos::RCP<const Tpetra_RowMatrix>& Matrix) :
  Matrix_(Matrix),
  NumRows_(0),
  NumNonzeros_(0),
  MaxNumEntries_(0),
  MaxNumEntriesA_(0)
{
  sprintf(Label_,"%s","Tifpack_LocalFilter");

#ifdef HAVE_MPI
  SerialComm_ = Teuchos::rcp( new Tpetra_MpiComm(MPI_COMM_SELF) );
#else
  SerialComm_ = Teuchos::rcp( new Tpetra_SerialComm );
#endif

  // localized matrix has all the local rows of Matrix
  NumRows_ = Matrix->NumMyRows();

  // build a linear map, based on the serial communicator
  Map_ = Teuchos::rcp( new Tpetra_Map(NumRows_,0,*SerialComm_) );

  // NumEntries_ will contain the actual number of nonzeros
  // for each localized row (that is, without external nodes,
  // and always with the diagonal entry)
  NumEntries_.resize(NumRows_);

  // want to store the diagonal vector. FIXME: am I really useful?
  Diagonal_ = Teuchos::rcp( new Tpetra_Vector(*Map_) );
  if (Diagonal_ == Teuchos::null) TIFPACK_CHK_ERRV(-5);
  
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
    TIFPACK_CHK_ERRV(ExtractMyRowCopy(i,MaxNumEntries_,Nnz,&Values_[0],&Indices_[0]));

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
int Tifpack_LocalFilter::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{
  if ((MyRow < 0) || (MyRow >= NumRows_)) {
    TIFPACK_CHK_ERR(-1); // range not valid
  }

  if (Length < NumEntries_[MyRow])
    return(-1);

  // always extract using the object Values_ and Indices_.
  // This is because I need more space than that given by
  // the user (for the external nodes)
  int Nnz;
  int ierr = Matrix_->ExtractMyRowCopy(MyRow,MaxNumEntriesA_,Nnz,
				       &Values_[0],&Indices_[0]);

  TIFPACK_CHK_ERR(ierr);

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
int Tifpack_LocalFilter::ExtractDiagonalCopy(Tpetra_Vector & Diagonal) const
{
  if (!Diagonal.Map().SameAs(*Map_))
    TIFPACK_CHK_ERR(-1);
  Diagonal = *Diagonal_;
  return(0);
}

//==============================================================================
int Tifpack_LocalFilter::Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
	  Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const 
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
    TIFPACK_CHK_ERR(ierr);

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
int Tifpack_LocalFilter::ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
		 Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  TIFPACK_CHK_ERR(-1); // not implemented
}

//==============================================================================
const Tpetra_BlockMap& Tifpack_LocalFilter::Map() const
{
  return(*Map_);
}
