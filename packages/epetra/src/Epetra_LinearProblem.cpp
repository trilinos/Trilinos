
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_LinearProblem.h"


//=============================================================================
Epetra_LinearProblem::Epetra_LinearProblem(void) 
  : A_(0),
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
  : A_(A),
    X_(X),
    B_(B),
    OperatorSymmetric_(false),
    PDL_(unsure),
    LeftScaled_(false),
    RightScaled_(false),
    LeftScaleVector_(0),
    RightScaleVector_(0)
{
}
//=============================================================================
Epetra_LinearProblem::Epetra_LinearProblem(const Epetra_LinearProblem& Problem) 
  : A_(Problem.A_),
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
  EPETRA_CHK_ERR(A_->LeftScale(D));
  EPETRA_CHK_ERR(B_->Multiply(1.0, D, *B_, 0.0));
  return(0);
}
//=============================================================================
int Epetra_LinearProblem::RightScale(const Epetra_Vector & D)
{
  EPETRA_CHK_ERR(A_->RightScale(D));
  EPETRA_CHK_ERR(X_->ReciprocalMultiply(1.0, D, *X_, 0.0));
  return(0);
}
//=============================================================================
int Epetra_LinearProblem::CheckInput() const {
  int ierr = 0;
  if (A_==0) ierr = -1;
  if (X_==0) ierr = -2;
  if (B_==0) ierr = -3;

  if (A_->NumGlobalRows()!=X_->GlobalLength()) ierr = -4;
  if (A_->NumGlobalCols()!=B_->GlobalLength()) ierr = -5;

  EPETRA_CHK_ERR(ierr);
  return(0);

}
