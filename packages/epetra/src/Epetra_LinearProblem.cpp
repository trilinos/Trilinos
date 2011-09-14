
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
