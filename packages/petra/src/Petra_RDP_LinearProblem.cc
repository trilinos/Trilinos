
#include "Petra_RDP_LinearProblem.h"


//=============================================================================
Petra_RDP_LinearProblem::Petra_RDP_LinearProblem(void) 
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
Petra_RDP_LinearProblem::Petra_RDP_LinearProblem(Petra_RDP_RowMatrix * A, 
					       Petra_RDP_MultiVector * X,
					       Petra_RDP_MultiVector * B) 
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
Petra_RDP_LinearProblem::Petra_RDP_LinearProblem(const Petra_RDP_LinearProblem& Problem) 
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
Petra_RDP_LinearProblem::~Petra_RDP_LinearProblem(void)  
{
}
//=============================================================================
int Petra_RDP_LinearProblem::LeftScale(const Petra_RDP_Vector & D)
{
  int ierr1 = 0, ierr2 = 0;
  ierr1 = A_->LeftScale(D);
  if (ierr1<0) PETRA_CHK_ERR(ierr1);
  ierr2 = B_->Multiply(1.0, D, *B_, 0.0);
  if (ierr2<0) PETRA_CHK_ERR(ierr2);
  
  
  if (ierr1>0) PETRA_CHK_ERR(ierr1);
  if (ierr2>0) PETRA_CHK_ERR(ierr2);

  return(0);
}
//=============================================================================
int Petra_RDP_LinearProblem::RightScale(const Petra_RDP_Vector & D)
{
  int ierr1 = 0, ierr2 = 0;
  ierr1 = A_->RightScale(D);
  if (ierr1<0) PETRA_CHK_ERR(ierr1);
  ierr2 = X_->ReciprocalMultiply(1.0, D, *X_, 0.0);
  if (ierr2<0) PETRA_CHK_ERR(ierr2);
  
  if (ierr1>0) PETRA_CHK_ERR(ierr1);
  if (ierr2>0) PETRA_CHK_ERR(ierr2);

  return(0);
}
//=============================================================================
int Petra_RDP_LinearProblem::CheckInput() const {
  int ierr = 0;
  if (A_==0) ierr = -1;
  if (X_==0) ierr = -2;
  if (B_==0) ierr = -3;

  if (A_->NumGlobalRows()!=X_->GlobalLength()) ierr = -4;
  if (A_->NumGlobalCols()!=B_->GlobalLength()) ierr = -5;

  PETRA_CHK_ERR(ierr);

}
