
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

#include "Epetra_SerialSpdDenseSolver.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_SerialDenseMatrix.h"

//=============================================================================
Epetra_SerialSpdDenseSolver::Epetra_SerialSpdDenseSolver(void)
  : Epetra_SerialDenseSolver(),
    SCOND_(-1.0),
    SymMatrix_(0),
    SymFactor_(0)
{
}
//=============================================================================
Epetra_SerialSpdDenseSolver::~Epetra_SerialSpdDenseSolver()
{
  if (SymFactor_ != SymMatrix_ && SymFactor_ != 0) {
    delete SymFactor_; SymFactor_ = 0; Factor_ = 0;
  }
}
//=============================================================================
int Epetra_SerialSpdDenseSolver::SetMatrix(Epetra_SerialSymDenseMatrix & A_in) {
  
  SymMatrix_=&A_in;
  SymFactor_=&A_in;
  SCOND_ = -1.0;
  // Also call SerialDensematrix set method
  return(Epetra_SerialDenseSolver::SetMatrix( (Epetra_SerialDenseMatrix &) A_in));
}
//=============================================================================
int Epetra_SerialSpdDenseSolver::Factor(void) {
  if (Factored()) return(0); // Already factored
  if (Inverted()) EPETRA_CHK_ERR(-100); // Cannot factor inverted matrix
  int ierr = 0;

  ANORM_ = SymMatrix_->OneNorm();


  // If we want to refine the solution, then the factor must
  // be stored separatedly from the original matrix

  if (A_ == AF_)
    if (RefineSolution_ ) {
      SymFactor_ = new Epetra_SerialSymDenseMatrix(*SymMatrix_);
      Factor_ = SymFactor_;
      AF_ = SymFactor_->A();
      LDAF_ = SymFactor_->LDA();
    }
  if (Equilibrate_) ierr = EquilibrateMatrix();

  if (ierr!=0) EPETRA_CHK_ERR(ierr-2);
  
  POTRF (SymMatrix_->UPLO(), N_, AF_, LDAF_, &INFO_);
  Factored_ = true;
  double DN = N_;
  UpdateFlops((DN*DN*DN)/3.0);

  EPETRA_CHK_ERR(INFO_);
  return(0);

}

//=============================================================================
int Epetra_SerialSpdDenseSolver::Solve(void) {
  int ierr = 0;

  // We will call one of four routines depending on what services the user wants and 
  // whether or not the matrix has been inverted or factored already.
  //
  // If the matrix has been inverted, use DGEMM to compute solution.
  // Otherwise, if the user want the matrix to be equilibrated or wants a refined solution, we will
  // call the X interface.
  // Otherwise, if the matrix is already factored we will call the TRS interface.
  // Otherwise, if the matrix is unfactored we will call the SV interface.

  if (Equilibrate_) {
    ierr = Epetra_SerialDenseSolver::EquilibrateRHS();
    B_Equilibrated_ = true;
  }
  EPETRA_CHK_ERR(ierr);
  if (A_Equilibrated_ && !B_Equilibrated_) EPETRA_CHK_ERR(-1); // Matrix and vectors must be similarly scaled
  if (!A_Equilibrated_ && B_Equilibrated_) EPETRA_CHK_ERR(-2);
  if (B_==0) EPETRA_CHK_ERR(-3); // No B
  if (X_==0) EPETRA_CHK_ERR(-4); // No B

  if (ShouldEquilibrate() && !A_Equilibrated_) ierr = 1; // Warn that the system should be equilibrated.

  double DN = N_;
  double DNRHS = NRHS_;
  if (Inverted()) {

    if (B_==X_) EPETRA_CHK_ERR(-100); // B and X must be different for this case
    GEMM('N', 'N', N_, NRHS_, N_, 1.0, AF_, LDAF_,
		B_, LDB_, 0.0, X_, LDX_);
    if (INFO_!=0) EPETRA_CHK_ERR(INFO_);
    UpdateFlops(2.0*DN*DN*DNRHS);
    Solved_ = true;
  }
  else {

    if (!Factored()) Factor(); // Matrix must be factored
    if (B_!=X_) {
       *LHS_ = *RHS_; // Copy B to X if needed 
       X_ = LHS_->A(); LDX_ = LHS_->LDA();
    }
    
    POTRS(SymMatrix_->UPLO(), N_, NRHS_, AF_, LDAF_, X_, LDX_, &INFO_);
    if (INFO_!=0) EPETRA_CHK_ERR(INFO_);
    UpdateFlops(2.0*DN*DN*DNRHS);
    Solved_ = true;

  }
  int ierr1=0;
  if (RefineSolution_) ierr1 = ApplyRefinement();
  if (ierr1!=0) {
    EPETRA_CHK_ERR(ierr1);
  }
  else {
    EPETRA_CHK_ERR(ierr);
  }
  if (Equilibrate_) ierr1 = Epetra_SerialDenseSolver::UnequilibrateLHS();
  EPETRA_CHK_ERR(ierr1);
  return(0);
}
//=============================================================================
int Epetra_SerialSpdDenseSolver::ApplyRefinement(void)
{
  double DN = N_;
  double DNRHS = NRHS_;
  if (!Solved()) EPETRA_CHK_ERR(-100); // Must have an existing solution
  if (A_==AF_) EPETRA_CHK_ERR(-101); // Cannot apply refine if no original copy of A.

  if (FERR_ != 0) delete [] FERR_; // Always start with a fresh copy of FERR_, since NRHS_ may change
  FERR_ = new double[NRHS_];
  if (BERR_ != 0) delete [] BERR_; // Always start with a fresh copy of FERR_, since NRHS_ may change
  BERR_ = new double[NRHS_];
  AllocateWORK();
  AllocateIWORK();
  
  PORFS(SymMatrix_->UPLO(), N_, NRHS_, A_, LDA_, AF_, LDAF_,
	       B_, LDB_, X_, LDX_, FERR_, BERR_, 
	       WORK_, IWORK_, &INFO_);
  
  
  SolutionErrorsEstimated_ = true;
  ReciprocalConditionEstimated_ = true;
  SolutionRefined_ = true;
  
  UpdateFlops(2.0*DN*DN*DNRHS); // Not sure of count
  
  EPETRA_CHK_ERR(INFO_);
  return(0);
  
}

//=============================================================================
int Epetra_SerialSpdDenseSolver::ComputeEquilibrateScaling(void) {
  if (R_!=0) return(0); // Already computed
 
  double DN = N_;
  R_ = new double[N_];
  C_ = R_;
  
  POEQU (N_, AF_, LDAF_, R_, &SCOND_, &AMAX_, &INFO_);
  if (INFO_ != 0) EPETRA_CHK_ERR(INFO_);

  if (SCOND_<0.1 || AMAX_ < Epetra_Underflow || AMAX_ > Epetra_Overflow) ShouldEquilibrate_ = true;

  C_ = R_; // Set column scaling pointer so we can use EquilibrateRHS and UnequilibrateLHS from base class
  UpdateFlops(2.0*DN*DN);
  
  return(0);
}

//=============================================================================
int Epetra_SerialSpdDenseSolver::EquilibrateMatrix(void) {
  int i, j;
  int ierr = 0;

  if (A_Equilibrated_) return(0); // Already done
  if (R_==0) ierr = ComputeEquilibrateScaling(); // Compute S if needed
  if (ierr!=0) EPETRA_CHK_ERR(ierr);
  if (SymMatrix_->Upper()) {
    if (A_==AF_) {
      double * ptr;
      for (j=0; j<N_; j++) {
	ptr = A_ + j*LDA_;
	double s1 = R_[j];
	for (i=0; i<=j; i++) {
	  *ptr = *ptr*s1*R_[i];
	  ptr++;
	}
      }
    }
    else {
      double * ptr;
      double * ptr1;
      for (j=0; j<N_; j++) {
	ptr = A_ + j*LDA_;
	ptr1 = AF_ + j*LDAF_;
	double s1 = R_[j];
	for (i=0; i<=j; i++) {
	  *ptr = *ptr*s1*R_[i];
	  ptr++;
	  *ptr1 = *ptr1*s1*R_[i];
	  ptr1++;
	}
      }
    }
  }
  else {
    if (A_==AF_) {
      double * ptr;
      for (j=0; j<N_; j++) {
	ptr = A_ + j + j*LDA_;
	double s1 = R_[j];
	for (i=j; i<N_; i++) {
	  *ptr = *ptr*s1*R_[i];
	  ptr++;
	}
      }
    }
    else {
      double * ptr;
      double * ptr1;
      for (j=0; j<N_; j++) {
	ptr = A_ + j + j*LDA_;
	ptr1 = AF_ + j + j*LDAF_;
	double s1 = R_[j];
	for (i=j; i<N_; i++) {
	  *ptr = *ptr*s1*R_[i];
	  ptr++;
	  *ptr1 = *ptr1*s1*R_[i];
	  ptr1++;
	}
      }
    }
  }
  A_Equilibrated_ = true;
  double NumFlops = (double) ((N_+1)*N_/2);
  if (A_==AF_) NumFlops += NumFlops;
  UpdateFlops(NumFlops);
  
  return(0);
}

//=============================================================================
int Epetra_SerialSpdDenseSolver::Invert(void)
{
  if (!Factored()) Factor(); // Need matrix factored.
  POTRI ( SymMatrix_->UPLO(), N_, AF_, LDAF_, &INFO_);
  // Copy lower/upper triangle to upper/lower triangle: make full inverse
  SymMatrix_->CopyUPLOMat(SymMatrix_->Upper(), AF_, LDAF_, N_);
  double DN = N_;
  UpdateFlops((DN*DN*DN));
  Inverted_ = true;
  Factored_ = false;
  
  EPETRA_CHK_ERR(INFO_);
  return(0);
}

//=============================================================================
int Epetra_SerialSpdDenseSolver::ReciprocalConditionEstimate(double & Value)
{
  int ierr = 0;
  if (ReciprocalConditionEstimated()) {
    Value = RCOND_;
    return(0); // Already computed, just return it.
  }

  if (ANORM_<0.0) ANORM_ = SymMatrix_->OneNorm();
  if (ierr!=0) EPETRA_CHK_ERR(ierr-1);
  if (!Factored()) ierr = Factor(); // Need matrix factored.
  if (ierr!=0) EPETRA_CHK_ERR(ierr-2);

  AllocateWORK();
  AllocateIWORK();
  POCON( SymMatrix_->UPLO(), N_, AF_, LDAF_, ANORM_, &RCOND_, WORK_, IWORK_, &INFO_);
  ReciprocalConditionEstimated_ = true;
  Value = RCOND_;
  UpdateFlops(2*N_*N_); // Not sure of count
  EPETRA_CHK_ERR(INFO_);
  return(0);
}
