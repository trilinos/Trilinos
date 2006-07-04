
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "Epetra_LinearProblem.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"


//=============================================================================
Epetra_LinearProblem::Epetra_LinearProblem(void) 
  : Operator_(0),
    A_(0),
    X_(0),
    B_(0),
    OperatorSymmetric_(false),
    PDL_(unsure),
    LeftScaled_(false),
    RightScaled_(false),
    LeftScaleVector_(0),
    RightScaleVector_(0)
{
}
//=============================================================================
Epetra_LinearProblem::Epetra_LinearProblem(Epetra_RowMatrix * A, 
					       Epetra_MultiVector * X,
					       Epetra_MultiVector * B) 
  : Operator_(0),
    A_(A),
    X_(X),
    B_(B),
    OperatorSymmetric_(false),
    PDL_(unsure),
    LeftScaled_(false),
    RightScaled_(false),
    LeftScaleVector_(0),
    RightScaleVector_(0)
{
  Operator_ = dynamic_cast<Epetra_Operator *>(A_); // Try to make matrix an operator
}
//=============================================================================
Epetra_LinearProblem::Epetra_LinearProblem(Epetra_Operator * A, 
					       Epetra_MultiVector * X,
					       Epetra_MultiVector * B) 
  : Operator_(A),
    A_(0),
    X_(X),
    B_(B),
    OperatorSymmetric_(false),
    PDL_(unsure),
    LeftScaled_(false),
    RightScaled_(false),
    LeftScaleVector_(0),
    RightScaleVector_(0)
{
  A_ = dynamic_cast<Epetra_RowMatrix *>(Operator_); // Try to make operator a matrix
}
//=============================================================================
Epetra_LinearProblem::Epetra_LinearProblem(const Epetra_LinearProblem& Problem) 
  : Operator_(Problem.Operator_),
    A_(Problem.A_),
    X_(Problem.X_),
    B_(Problem.B_),
    OperatorSymmetric_(Problem.OperatorSymmetric_),
    PDL_(Problem.PDL_),
    LeftScaled_(Problem.LeftScaled_),
    RightScaled_(Problem.RightScaled_),
    LeftScaleVector_(Problem.LeftScaleVector_),
    RightScaleVector_(Problem.RightScaleVector_)
{
}
//=============================================================================
Epetra_LinearProblem::~Epetra_LinearProblem(void)  
{
}
//=============================================================================
int Epetra_LinearProblem::LeftScale(const Epetra_Vector & D)
{
  if (A_==0) EPETRA_CHK_ERR(-1); // No matrix defined
  if (B_==0) EPETRA_CHK_ERR(-2); // No RHS defined
  if (A_->UseTranspose()) {
    EPETRA_CHK_ERR(A_->RightScale(D));
    EPETRA_CHK_ERR(X_->ReciprocalMultiply(1.0, D, *X_, 0.0));
  }
  else {
    EPETRA_CHK_ERR(A_->LeftScale(D));
    EPETRA_CHK_ERR(B_->Multiply(1.0, D, *B_, 0.0));
  }

  return(0);
}
//=============================================================================
int Epetra_LinearProblem::RightScale(const Epetra_Vector & D)
{
  if (A_==0) EPETRA_CHK_ERR(-1); // No matrix defined
  if (X_==0) EPETRA_CHK_ERR(-2); // No LHS defined
  if (A_->UseTranspose()) {
    EPETRA_CHK_ERR(A_->LeftScale(D));
    EPETRA_CHK_ERR(B_->Multiply(1.0, D, *B_, 0.0));
  }
  else {
    EPETRA_CHK_ERR(A_->RightScale(D));
    EPETRA_CHK_ERR(X_->ReciprocalMultiply(1.0, D, *X_, 0.0));
  }
  return(0);
}
//=============================================================================
int Epetra_LinearProblem::CheckInput() const {
  int ierr = 0;
  if (Operator_==0) ierr = -1;
  if (X_==0) ierr = -2;
  if (B_==0) ierr = -3;

  EPETRA_CHK_ERR(ierr);  // Return now if any essential objects missing

  if (A_==0) EPETRA_CHK_ERR(1); // Return warning error because this problem has no matrix (just an operator)

  if (!A_->OperatorDomainMap().SameAs(X_->Map())) ierr = -4;
  if (!A_->OperatorRangeMap().SameAs(B_->Map())) ierr = -5;

  EPETRA_CHK_ERR(ierr);
  return(0);

}
