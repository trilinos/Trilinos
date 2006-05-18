
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

#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseMatrix.h"

//=============================================================================
Epetra_SerialDenseSolver::Epetra_SerialDenseSolver()
  : Epetra_CompObject(),
    Epetra_BLAS(),
    Epetra_LAPACK(),
    Equilibrate_(false),
    ShouldEquilibrate_(false),
    A_Equilibrated_(false),
    B_Equilibrated_(false),
    Transpose_(false),
    Factored_(false),
    EstimateSolutionErrors_(false),
    SolutionErrorsEstimated_(false),
    Solved_(false),
    Inverted_(false),
    ReciprocalConditionEstimated_(false),
    RefineSolution_(false),
    SolutionRefined_(false),
    TRANS_('N'),
    M_(0),
    N_(0),
    Min_MN_(0),
    NRHS_(0),
    LDA_(0),
    LDAF_(0),
    LDB_(0),
    LDX_(0),
    INFO_(0),
    LWORK_(0),
    IPIV_(0),
    IWORK_(0),
    ANORM_(0.0),
    RCOND_(0.0),
    ROWCND_(0.0),
    COLCND_(0.0),
    AMAX_(0.0),
    Matrix_(0),
    LHS_(0),
    RHS_(0),
    Factor_(0),
    A_(0),
    FERR_(0),
    BERR_(0),
    AF_(0),
    WORK_(0),
    R_(0),
    C_(0),
    B_(0),
    X_(0)
{
  InitPointers();
  ResetMatrix();
  ResetVectors();
}
//=============================================================================
Epetra_SerialDenseSolver::~Epetra_SerialDenseSolver()
{
  DeleteArrays();
}
//=============================================================================
void Epetra_SerialDenseSolver::InitPointers()
{
  IWORK_ = 0;
  FERR_ = 0;
  BERR_ = 0;
  Factor_ =0;
  Matrix_ =0;
  AF_ = 0;
  IPIV_ = 0;
  WORK_ = 0;
  R_ = 0;
  C_ = 0;
  INFO_ = 0;
  LWORK_ = 0;    
}
//=============================================================================
void Epetra_SerialDenseSolver::DeleteArrays()
{
  if (IWORK_ != 0) {delete [] IWORK_; IWORK_ = 0;}
  if (FERR_ != 0)  {delete [] FERR_; FERR_ = 0;}
  if (BERR_ != 0)  {delete [] BERR_; BERR_ = 0;}
  if (Factor_ != Matrix_ && Factor_ != 0)   {delete Factor_; Factor_ = 0;}
  if (Factor_ !=0) Factor_ = 0;
  if (AF_ !=0) AF_ = 0;
  if (IPIV_ != 0)  {delete [] IPIV_;IPIV_ = 0;}
  if (WORK_ != 0)  {delete [] WORK_;WORK_ = 0;}
  if (R_ != 0 && R_ != C_)     {delete [] R_;R_ = 0;}
  if (R_ != 0) R_ = 0;
  if (C_ != 0)     {delete [] C_;C_ = 0;}
  INFO_ = 0;
  LWORK_ = 0;    
}
//=============================================================================
void Epetra_SerialDenseSolver::ResetMatrix()
{
  DeleteArrays();
  ResetVectors();
  Matrix_ = 0;
  Factor_ = 0;
  A_Equilibrated_ = false;
  Factored_ = false;
  Inverted_ = false;
  M_ = 0;
  N_ = 0;
  Min_MN_ = 0;
  LDA_ = 0;
  LDAF_ = 0;
  ANORM_ = -1.0;
  RCOND_ = -1.0;
  ROWCND_ = -1.0;
  COLCND_ = -1.0;
  AMAX_ = -1.0;
  A_ = 0;

}
//=============================================================================
int Epetra_SerialDenseSolver::SetMatrix(Epetra_SerialDenseMatrix & A) {
  ResetMatrix();
  Matrix_ = &A;
  Factor_ = &A;
  M_ = A.M();
  N_ = A.N();
  Min_MN_ = EPETRA_MIN(M_,N_);
  LDA_ = A.LDA();
  LDAF_ = LDA_;
  A_ = A.A();
  AF_ = A.A();
  return(0);
}
//=============================================================================
void Epetra_SerialDenseSolver::ResetVectors()
{
  LHS_ = 0;
  RHS_ = 0;
  B_ = 0;
  X_ = 0;
  ReciprocalConditionEstimated_ = false;
  SolutionRefined_ = false;
  Solved_ = false;
  SolutionErrorsEstimated_ = false;
  B_Equilibrated_ = false;
  NRHS_ = 0;
  LDB_ = 0;
  LDX_ = 0;
}
//=============================================================================
int Epetra_SerialDenseSolver::SetVectors(Epetra_SerialDenseMatrix & X, Epetra_SerialDenseMatrix & B)
{
  if (B.M()!=X.M() || B.N() != X.N()) EPETRA_CHK_ERR(-1);
  if (B.A()==0) EPETRA_CHK_ERR(-2);
  if (B.LDA()<1) EPETRA_CHK_ERR(-3);
  if (X.A()==0) EPETRA_CHK_ERR(-4);
  if (X.LDA()<1) EPETRA_CHK_ERR(-5);

  ResetVectors(); 
  LHS_ = &X;
  RHS_ = &B;
  NRHS_ = B.N();

  B_ = B.A();
  LDB_ = B.LDA();
  X_ = X.A();
  LDX_ = X.LDA();
  return(0);
}
//=============================================================================
void Epetra_SerialDenseSolver::EstimateSolutionErrors(bool Flag) {
  EstimateSolutionErrors_ = Flag;
  // If the errors are estimated, this implies that the solution must be refined
  RefineSolution_ = RefineSolution_ || Flag;
  return;
}
//=============================================================================
int Epetra_SerialDenseSolver::Factor(void) {
  if (Factored()) return(0); // Already factored
  if (Inverted()) EPETRA_CHK_ERR(-100); // Cannot factor inverted matrix
  int ierr = 0;

  ANORM_ = Matrix_->OneNorm(); // Compute 1-Norm of A


  // If we want to refine the solution, then the factor must
  // be stored separatedly from the original matrix

  if (A_ == AF_)
    if (RefineSolution_ ) {
      Factor_ = new Epetra_SerialDenseMatrix(*Matrix_);
      AF_ = Factor_->A();
      LDAF_ = Factor_->LDA();
    }
  
  if (Equilibrate_) ierr = EquilibrateMatrix();

  if (ierr!=0) EPETRA_CHK_ERR(ierr-2);
  
  if (IPIV_==0) IPIV_ = new int[Min_MN_]; // Allocated Pivot vector if not already done.

  GETRF (M_, N_, AF_, LDAF_, IPIV_, &INFO_);

  Factored_ = true;
  double DN = N_;
  UpdateFlops(2.0*(DN*DN*DN)/3.0);

  EPETRA_CHK_ERR(INFO_);
  return(0);

}

//=============================================================================
int Epetra_SerialDenseSolver::Solve(void) {
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
    ierr = EquilibrateRHS();
    B_Equilibrated_ = true;
  }
  EPETRA_CHK_ERR(ierr);
  if (A_Equilibrated_ && !B_Equilibrated_) EPETRA_CHK_ERR(-1); // Matrix and vectors must be similarly scaled
  if (!A_Equilibrated_ && B_Equilibrated_) EPETRA_CHK_ERR(-2);
  if (B_==0) EPETRA_CHK_ERR(-3); // No B
  if (X_==0) EPETRA_CHK_ERR(-4); // No X

  if (ShouldEquilibrate() && !A_Equilibrated_) ierr = 1; // Warn that the system should be equilibrated.

  double DN = N_;
  double DNRHS = NRHS_;
  if (Inverted()) {

    if (B_==X_) EPETRA_CHK_ERR(-100); // B and X must be different for this case

    GEMM(TRANS_, 'N', N_, NRHS_, N_, 1.0, AF_, LDAF_,
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
    GETRS(TRANS_, N_, NRHS_, AF_, LDAF_, IPIV_, X_, LDX_, &INFO_);
    if (INFO_!=0) EPETRA_CHK_ERR(INFO_);
    UpdateFlops(2.0*DN*DN*DNRHS);
    Solved_ = true;

  }
  int ierr1=0;
  if (RefineSolution_ && !Inverted()) ierr1 = ApplyRefinement();
  if (ierr1!=0) EPETRA_CHK_ERR(ierr1)
  else
    EPETRA_CHK_ERR(ierr);
  
  if (Equilibrate_) ierr1 = UnequilibrateLHS();
  EPETRA_CHK_ERR(ierr1);
  return(0);
}
//=============================================================================
int Epetra_SerialDenseSolver::ApplyRefinement(void)
{
  double DN = N_;
  double DNRHS = NRHS_;
  if (!Solved()) EPETRA_CHK_ERR(-100); // Must have an existing solution
  if (A_==AF_) EPETRA_CHK_ERR(-101); // Cannot apply refine if no original copy of A.

  if (FERR_ != 0) delete [] FERR_; // Always start with a fresh copy of FERR_, since NRHS_ may change
  FERR_ = new double[NRHS_];
  if (BERR_ != 0) delete [] BERR_; // Always start with a fresh copy of BERR_, since NRHS_ may change
  BERR_ = new double[NRHS_];
  AllocateWORK();
  AllocateIWORK();
  
  GERFS(TRANS_, N_, NRHS_, A_, LDA_, AF_, LDAF_, IPIV_,
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
int Epetra_SerialDenseSolver::ComputeEquilibrateScaling(void) {
  if (R_!=0) return(0); // Already computed
 
  double DM = M_;
  double DN = N_;
  R_ = new double[M_];
  C_ = new double[N_];
  
  GEEQU (M_, N_, AF_, LDAF_, R_, C_, &ROWCND_, &COLCND_, &AMAX_, &INFO_);
  if (INFO_ != 0) EPETRA_CHK_ERR(INFO_);

  if (COLCND_<0.1 || ROWCND_<0.1 || AMAX_ < Epetra_Underflow || AMAX_ > Epetra_Overflow) ShouldEquilibrate_ = true;

  UpdateFlops(4.0*DM*DN);
  
  return(0);
}

//=============================================================================
int Epetra_SerialDenseSolver::EquilibrateMatrix(void)
{
  int i, j;
  int ierr = 0;

  double DN = N_;
  double DM = M_;

  if (A_Equilibrated_) return(0); // Already done
  if (R_==0) ierr = ComputeEquilibrateScaling(); // Compute R and C if needed
  if (ierr!=0) EPETRA_CHK_ERR(ierr);
  if (A_==AF_) {
    double * ptr;
    for (j=0; j<N_; j++) {
      ptr = A_ + j*LDA_;
      double s1 = C_[j];
      for (i=0; i<M_; i++) {
	*ptr = *ptr*s1*R_[i];
	ptr++;
      }
    }
    UpdateFlops(2.0*DM*DN);
  }
  else {
    double * ptr;
    double * ptr1;
    for (j=0; j<N_; j++) {
      ptr = A_ + j*LDA_;
      ptr1 = AF_ + j*LDAF_;
      double s1 = C_[j];
      for (i=0; i<M_; i++) {
	*ptr = *ptr*s1*R_[i];
	ptr++;
	*ptr1 = *ptr1*s1*R_[i];
	ptr1++;
      }
    }
    UpdateFlops(4.0*DM*DN);
  }
  
  A_Equilibrated_ = true;
  
  return(0);
}

//=============================================================================
int Epetra_SerialDenseSolver::EquilibrateRHS(void)
{
  int i, j;
  int ierr = 0;

  if (B_Equilibrated_) return(0); // Already done
  if (R_==0) ierr = ComputeEquilibrateScaling(); // Compute R and C if needed
  if (ierr!=0) EPETRA_CHK_ERR(ierr);

  double * R = R_;
  if (Transpose_) R = C_;

  double * ptr;
  for (j=0; j<NRHS_; j++) {
    ptr = B_ + j*LDB_;
    for (i=0; i<M_; i++) {
      *ptr = *ptr*R[i];
      ptr++;
    }
  }

  
  B_Equilibrated_ = true;
  UpdateFlops((double) N_*(double) NRHS_);
  
  return(0);
}

//=============================================================================
int Epetra_SerialDenseSolver::UnequilibrateLHS(void)
{
  int i, j;

  if (!B_Equilibrated_) return(0); // Nothing to do

  double * C = C_;
  if (Transpose_) C = R_;

  double * ptr;
  for (j=0; j<NRHS_; j++) {
    ptr = X_ + j*LDX_;
    for (i=0; i<N_; i++) {
      *ptr = *ptr*C[i];
      ptr++;
    }
  }

  
  UpdateFlops((double) N_ *(double) NRHS_);
  
  return(0);
}

//=============================================================================
int Epetra_SerialDenseSolver::Invert(void)
{
  if (!Factored()) Factor(); // Need matrix factored.

  /* This section work with LAPACK Version 3.0 only 
  // Setting LWORK = -1 and calling GETRI will return optimal work space size in WORK_TMP
  int LWORK_TMP = -1;
  double WORK_TMP;
  GETRI ( N_, AF_, LDAF_, IPIV_, &WORK_TMP, &LWORK_TMP, &INFO_);
  LWORK_TMP = WORK_TMP; // Convert to integer
  if (LWORK_TMP>LWORK_) {
  if (WORK_!=0) delete [] WORK_;
  LWORK_ = LWORK_TMP;
  WORK_ = new double[LWORK_];
  }
  */
  // This section will work with any version of LAPACK 
  AllocateWORK();

  GETRI ( N_, AF_, LDAF_, IPIV_, WORK_, &LWORK_, &INFO_);

  double DN = N_;
  UpdateFlops((DN*DN*DN));
  Inverted_ = true;
  Factored_ = false;
  
  EPETRA_CHK_ERR(INFO_);
  return(0);
}

//=============================================================================
int Epetra_SerialDenseSolver::ReciprocalConditionEstimate(double & Value)
{
  int ierr = 0;
  if (ReciprocalConditionEstimated()) {
    Value = RCOND_;
    return(0); // Already computed, just return it.
  }

  if (ANORM_<0.0) ANORM_ = Matrix_->OneNorm();
  if (!Factored()) ierr = Factor(); // Need matrix factored.
  if (ierr!=0) EPETRA_CHK_ERR(ierr-2);

  AllocateWORK();
  AllocateIWORK();
  // We will assume a one-norm condition number
  GECON( '1', N_, AF_, LDAF_, ANORM_, &RCOND_, WORK_, IWORK_, &INFO_);
  ReciprocalConditionEstimated_ = true;
  Value = RCOND_;
  UpdateFlops(2*N_*N_); // Not sure of count
  EPETRA_CHK_ERR(INFO_);
  return(0);
}
//=============================================================================
void Epetra_SerialDenseSolver::Print(ostream& os) const {

  if (Matrix_!=0) os << *Matrix_;
  if (Factor_!=0) os << *Factor_;
  if (LHS_!=0) os << *LHS_;
  if (RHS_!=0) os << *RHS_;

}
