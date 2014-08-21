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
#include "Ifpack_TriDiContainer.h"
#include "Epetra_RowMatrix.h"

//==============================================================================
int Ifpack_TriDiContainer::NumRows() const
{
  return(NumRows_);
}

//==============================================================================
int Ifpack_TriDiContainer::Initialize()
{
  
  IsInitialized_ = false;

  IFPACK_CHK_ERR(LHS_.Reshape(NumRows_,NumVectors_));
  IFPACK_CHK_ERR(RHS_.Reshape(NumRows_,NumVectors_));
  IFPACK_CHK_ERR(ID_.Reshape(NumRows_,NumVectors_));
  IFPACK_CHK_ERR(Matrix_.Reshape(NumRows_,NumRows_));

  // zero out matrix elements
  // for (int i = 0 ; i < NumRows_ ; ++i)
  //   for (int j = 0 ; j < NumRows_ ; ++j)
  //     Matrix_(i,j) = 0.0;
  memset(Matrix_.A(),0,sizeof(double)*4*(NumRows_-1));

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
double& Ifpack_TriDiContainer::LHS(const int i, const int Vector)
{
  return(LHS_.A()[Vector * NumRows_ + i]);
}
  
//==============================================================================
double& Ifpack_TriDiContainer::RHS(const int i, const int Vector)
{
  return(RHS_.A()[Vector * NumRows_ + i]);
}

//==============================================================================
int Ifpack_TriDiContainer::
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
int Ifpack_TriDiContainer::ApplyInverse()
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
int& Ifpack_TriDiContainer::ID(const int i)
{
  return(ID_[i]);
}

//==============================================================================
// FIXME: optimize performances of this guy...
int Ifpack_TriDiContainer::Extract(const Epetra_RowMatrix& Matrix_in)
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
  std::vector<double> Values;
  Values.resize(Length);
  std::vector<int> Indices;
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
int Ifpack_TriDiContainer::Compute(const Epetra_RowMatrix& Matrix_in)
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

  Label_ = "Ifpack_TriDiContainer";

  // not sure of count
#ifdef IFPACK_FLOPCOUNTERS
  ComputeFlops_ += 4.0 * NumRows_ * NumRows_ * NumRows_ / 3;
#endif
  IsComputed_ = true;

  return(0);
}

//==============================================================================
 int Ifpack_TriDiContainer::Apply()
 {
   IFPACK_CHK_ERR(-300);
   //   if (IsComputed() == false)
   //     IFPACK_CHK_ERR(-3);

   //   if (KeepNonFactoredMatrix_) {
   //     IFPACK_CHK_ERR(RHS_.Multiply('N','N', 1.0,NonFactoredMatrix_,LHS_,0.0));
   //   }
   //   else
   //     IFPACK_CHK_ERR(RHS_.Multiply('N','N', 1.0,Matrix_,LHS_,0.0));

   // #ifdef IFPACK_FLOPCOUNTERS
   //   ApplyFlops_ += 2 * NumRows_ * NumRows_;
   // #endif
   return(0);
 }

//==============================================================================
ostream& Ifpack_TriDiContainer::Print(ostream & os) const
{
    os << "================================================================================" << endl;
  os << "Ifpack_TriDiContainer" << endl;
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
