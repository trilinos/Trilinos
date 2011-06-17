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
#include "Ifpack_DenseContainer.h"
#include "Epetra_RowMatrix.h"

//==============================================================================
int Ifpack_DenseContainer::NumRows() const
{
  return(NumRows_);
}

//==============================================================================
int Ifpack_DenseContainer::Initialize()
{
  
  IsInitialized_ = false;

  IFPACK_CHK_ERR(LHS_.Reshape(NumRows_,NumVectors_));
  IFPACK_CHK_ERR(RHS_.Reshape(NumRows_,NumVectors_));
  IFPACK_CHK_ERR(ID_.Reshape(NumRows_,NumVectors_));
  IFPACK_CHK_ERR(Matrix_.Reshape(NumRows_,NumRows_));

  // zero out matrix elements
  for (int i = 0 ; i < NumRows_ ; ++i)
    for (int j = 0 ; j < NumRows_ ; ++j)
      Matrix_(i,j) = 0.0;

  // zero out vector elements
  for (int i = 0 ; i < NumRows_ ; ++i)
    for (int j = 0 ; j < NumVectors_ ; ++j) {
      LHS_(i,j) = 0.0;
      RHS_(i,j) = 0.0;
    }

  // Set to -1 ID_'s
  for (int i = 0 ; i < NumRows_ ; ++i)
    ID_(i) = -1;  

  if (NumRows_ != 0) {
    IFPACK_CHK_ERR(Solver_.SetMatrix(Matrix_));
    IFPACK_CHK_ERR(Solver_.SetVectors(LHS_,RHS_));
  }

  IsInitialized_ = true;
  return(0);
  
}

//==============================================================================
double& Ifpack_DenseContainer::LHS(const int i, const int Vector)
{
  return(LHS_.A()[Vector * NumRows_ + i]);
}
  
//==============================================================================
double& Ifpack_DenseContainer::RHS(const int i, const int Vector)
{
  return(RHS_.A()[Vector * NumRows_ + i]);
}

//==============================================================================
int Ifpack_DenseContainer::
SetMatrixElement(const int row, const int col, const double value)
{
  if (IsInitialized() == false)
    IFPACK_CHK_ERR(Initialize());

  if ((row < 0) || (row >= NumRows())) {
    IFPACK_CHK_ERR(-2); // not in range
  }

  if ((col < 0) || (col >= NumRows())) {
    IFPACK_CHK_ERR(-2); // not in range
  }

  Matrix_(row, col) = value;

  return(0);

}

//==============================================================================
int Ifpack_DenseContainer::ApplyInverse()
{

  if (!IsComputed()) {
    IFPACK_CHK_ERR(-1);
  }
  
  if (NumRows_ != 0)
    IFPACK_CHK_ERR(Solver_.Solve());

#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += 2.0 * NumVectors_ * NumRows_ * NumRows_;
#endif
  return(0);
}

//==============================================================================
int& Ifpack_DenseContainer::ID(const int i)
{
  return(ID_[i]);
}

//==============================================================================
// FIXME: optimize performances of this guy...
int Ifpack_DenseContainer::Extract(const Epetra_RowMatrix& Matrix_in)
{

  for (int j = 0 ; j < NumRows_ ; ++j) {
    // be sure that the user has set all the ID's
    if (ID(j) == -1)
      IFPACK_CHK_ERR(-2);
    // be sure that all are local indices
    if (ID(j) > Matrix_in.NumMyRows())
      IFPACK_CHK_ERR(-2);
  }

  // allocate storage to extract matrix rows.
  int Length = Matrix_in.MaxNumEntries();
  vector<double> Values;
  Values.resize(Length);
  vector<int> Indices;
  Indices.resize(Length);

  for (int j = 0 ; j < NumRows_ ; ++j) {

    int LRID = ID(j);

    int NumEntries;

    int ierr = 
      Matrix_in.ExtractMyRowCopy(LRID, Length, NumEntries, 
			      &Values[0], &Indices[0]);
    IFPACK_CHK_ERR(ierr);

    for (int k = 0 ; k < NumEntries ; ++k) {

      int LCID = Indices[k];

      // skip off-processor elements
      if (LCID >= Matrix_in.NumMyRows()) 
	continue;

      // for local column IDs, look for each ID in the list
      // of columns hosted by this object
      // FIXME: use STL
      int jj = -1;
      for (int kk = 0 ; kk < NumRows_ ; ++kk)
	if (ID(kk) == LCID)
	  jj = kk;

      if (jj != -1)
	SetMatrixElement(j,jj,Values[k]);

    }
  }

  return(0);
}

//==============================================================================
int Ifpack_DenseContainer::Compute(const Epetra_RowMatrix& Matrix_in)
{
  IsComputed_ = false;
  if (IsInitialized() == false) {
    IFPACK_CHK_ERR(Initialize());
  }

  if (KeepNonFactoredMatrix_)
    NonFactoredMatrix_ = Matrix_;

  // extract local rows and columns
  IFPACK_CHK_ERR(Extract(Matrix_in));

  if (KeepNonFactoredMatrix_)
    NonFactoredMatrix_ = Matrix_;

  // factorize the matrix using LAPACK
  if (NumRows_ != 0)
    IFPACK_CHK_ERR(Solver_.Factor());

  Label_ = "Ifpack_DenseContainer";

  // not sure of count
#ifdef IFPACK_FLOPCOUNTERS
  ComputeFlops_ += 4.0 * NumRows_ * NumRows_ * NumRows_ / 3;
#endif
  IsComputed_ = true;

  return(0);
}

//==============================================================================
int Ifpack_DenseContainer::Apply()
{
  if (IsComputed() == false)
    IFPACK_CHK_ERR(-3);

  if (KeepNonFactoredMatrix_) {
    IFPACK_CHK_ERR(RHS_.Multiply('N','N', 1.0,NonFactoredMatrix_,LHS_,0.0));
  }
  else
    IFPACK_CHK_ERR(RHS_.Multiply('N','N', 1.0,Matrix_,LHS_,0.0));

#ifdef IFPACK_FLOPCOUNTERS
  ApplyFlops_ += 2 * NumRows_ * NumRows_;
#endif
  return(0);
}

//==============================================================================
ostream& Ifpack_DenseContainer::Print(ostream & os) const
{
    os << "================================================================================" << endl;
  os << "Ifpack_DenseContainer" << endl;
  os << "Number of rows          = " << NumRows() << endl;
  os << "Number of vectors       = " << NumVectors() << endl;
  os << "IsInitialized()         = " << IsInitialized() << endl;
  os << "IsComputed()            = " << IsComputed() << endl;
#ifdef IFPACK_FLOPCOUNTERS
  os << "Flops in Compute()      = " << ComputeFlops() << endl; 
  os << "Flops in ApplyInverse() = " << ApplyInverseFlops() << endl; 
#endif
  os << "================================================================================" << endl;
  os << endl;

  return(os);
}
