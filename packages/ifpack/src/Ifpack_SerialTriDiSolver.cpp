
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

#include "Ifpack_SerialTriDiSolver.h"
#include "Ifpack_SerialTriDiMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include <iostream>

//=============================================================================
Ifpack_SerialTriDiSolver::Ifpack_SerialTriDiSolver()
  : Epetra_CompObject(),
    Epetra_BLAS(),
    // Equilibrate_(false),
    // ShouldEquilibrate_(false),
    // A_Equilibrated_(false),
    // B_Equilibrated_(false),
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
    N_(0),
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
    // R_(0),
    // C_(0),
    B_(0),
    X_(0)
{
  InitPointers();
  ResetMatrix();
  ResetVectors();
}
//=============================================================================
Ifpack_SerialTriDiSolver::~Ifpack_SerialTriDiSolver()
{
  DeleteArrays();
}
//=============================================================================
void Ifpack_SerialTriDiSolver::InitPointers()
{
  IWORK_ = 0;
  FERR_ = 0;
  BERR_ = 0;
  Factor_ =0;
  Matrix_ =0;
  AF_ = 0;
  IPIV_ = 0;
  WORK_ = 0;
  // R_ = 0;
  // C_ = 0;
  INFO_ = 0;
  LWORK_ = 0;
}
//=============================================================================
void Ifpack_SerialTriDiSolver::DeleteArrays()
{
  if (IWORK_ != 0) {delete [] IWORK_; IWORK_ = 0;}
  if (FERR_ != 0)  {delete [] FERR_; FERR_ = 0;}
  if (BERR_ != 0)  {delete [] BERR_; BERR_ = 0;}
  if (Factor_ != Matrix_ && Factor_ != 0)   {delete Factor_; Factor_ = 0;}
  if (Factor_ !=0) Factor_ = 0;
  if (AF_ !=0) AF_ = 0;
  if (IPIV_ != 0)  {delete [] IPIV_;IPIV_ = 0;}
  if (WORK_ != 0)  {delete [] WORK_;WORK_ = 0;}
  // if (R_ != 0 && R_ != C_)     {delete [] R_;R_ = 0;}
  // if (R_ != 0) R_ = 0;
  // if (C_ != 0)     {delete [] C_;C_ = 0;}
  INFO_ = 0;
  LWORK_ = 0;
}
//=============================================================================
void Ifpack_SerialTriDiSolver::ResetMatrix()
{
  DeleteArrays();
  ResetVectors();
  Matrix_ = 0;
  Factor_ = 0;
  // A_Equilibrated_ = false;
  Factored_ = false;
  Inverted_ = false;
  N_ = 0;
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
int Ifpack_SerialTriDiSolver::SetMatrix(Ifpack_SerialTriDiMatrix & A_in) {
  ResetMatrix();
  Matrix_ = &A_in;
  Factor_ = &A_in;
  N_ = A_in.N();
  A_ = A_in.A();
  LDA_ = A_in.LDA();
  LDAF_ = LDA_;
  AF_ = A_in.A();
  return(0);
}
//=============================================================================
void Ifpack_SerialTriDiSolver::ResetVectors()
{
  LHS_ = 0;
  RHS_ = 0;
  B_ = 0;
  X_ = 0;
  ReciprocalConditionEstimated_ = false;
  SolutionRefined_ = false;
  Solved_ = false;
  SolutionErrorsEstimated_ = false;
  // B_Equilibrated_ = false;
  NRHS_ = 0;
  LDB_ = 0;
  LDX_ = 0;
}
//=============================================================================
int Ifpack_SerialTriDiSolver::SetVectors(Epetra_SerialDenseMatrix & X_in, Epetra_SerialDenseMatrix & B_in)
{
  if (B_in.N() != X_in.N()) EPETRA_CHK_ERR(-1);
  if (B_in.A()==0) EPETRA_CHK_ERR(-2);
  if (X_in.A()==0) EPETRA_CHK_ERR(-4);

  ResetVectors();
  LHS_ = &X_in;
  RHS_ = &B_in;
  NRHS_ = B_in.N();

  B_ = B_in.A();
  X_ = X_in.A();
  return(0);
}
//=============================================================================
void Ifpack_SerialTriDiSolver::EstimateSolutionErrors(bool Flag) {
  EstimateSolutionErrors_ = Flag;
  // If the errors are estimated, this implies that the solution must be refined
  RefineSolution_ = RefineSolution_ || Flag;
  return;
}
//=============================================================================
int Ifpack_SerialTriDiSolver::Factor(void) {
  if (Factored()) return(0); // Already factored
  if (Inverted()) EPETRA_CHK_ERR(-100); // Cannot factor inverted matrix

  ANORM_ = Matrix_->OneNorm(); // Compute 1-Norm of A

  // If we want to refine the solution, then the factor must
  // be stored separatedly from the original matrix

  Ifpack_SerialTriDiMatrix * F = Matrix_;

  if (A_ == AF_)
    if (RefineSolution_ ) {
      Factor_ = new Ifpack_SerialTriDiMatrix(*Matrix_);
      F = Factor_;
      AF_ = Factor_->A();
      LDAF_ = Factor_->LDA();
    }

  if (IPIV_==0) IPIV_ = new int[N_]; // Allocated Pivot vector if not already done.
  
  double * DL_  = F->DL();
  double * D_   = F->D();
  double * DU_  = F->DU();
  double * DU2_ = F->DU2();

  lapack.GTTRF(N_, DL_, D_, DU_, DU2_, IPIV_, &INFO_);

  Factored_ = true;
  double DN = N_;
  UpdateFlops(2.0*(DN*DN*DN)/3.0);   // \todo This is very likely wrong...

  EPETRA_CHK_ERR(INFO_);
  return(0);

}

//=============================================================================
int Ifpack_SerialTriDiSolver::Solve(void) {
  int ierr = 0;

  // We will call one of four routines depending on what services the user wants and
  // whether or not the matrix has been inverted or factored already.
  //
  // If the matrix has been inverted, use DGEMM to compute solution.
  // Otherwise, if the user want the matrix to be equilibrated or wants a refined solution, we will
  // call the X interface.
  // Otherwise, if the matrix is already factored we will call the TRS interface.
  // Otherwise, if the matrix is unfactored we will call the SV interface.

  if (B_==0) EPETRA_CHK_ERR(-3); // No B
  if (X_==0) EPETRA_CHK_ERR(-4); // No X

  double DN = N_;
  double DNRHS = NRHS_;
  if (Inverted()) {

    EPETRA_CHK_ERR(-101);  // don't allow this \cbl

  }
  else {

    if (!Factored()) Factor(); // Matrix must be factored

    if (B_!=X_) {
       *LHS_ = *RHS_; // Copy B to X if needed
       X_ = LHS_->A(); 
    }

    Ifpack_SerialTriDiMatrix * F;
    if(A_ == AF_)
      F = Matrix_;
    else
      F = Factor_;

    double * DL_  = F->DL();
    double * D_   = F->D();
    double * DU_  = F->DU();
    double * DU2_ = F->DU2();

    lapack.GTTRS(TRANS_,N_,NRHS_,DL_,D_,DU_,DU2_,IPIV_,X_,N_,&INFO_);

    if (INFO_!=0) EPETRA_CHK_ERR(INFO_);
    UpdateFlops(2.0*DN*DN*DNRHS);
    Solved_ = true;

  }
  int ierr1=0;
  if (RefineSolution_ && !Inverted()) ierr1 = ApplyRefinement();
  if (ierr1!=0) EPETRA_CHK_ERR(ierr1)
  else
    EPETRA_CHK_ERR(ierr);

  // if (Equilibrate_) ierr1 = UnequilibrateLHS();
  // EPETRA_CHK_ERR(ierr1);
  return(0);
}
//=============================================================================
 int Ifpack_SerialTriDiSolver::ApplyRefinement(void)
 {
   std::cout<<" SerialTriDiSolver::ApplyRefinement this function is not supported"<<std::endl;
   EPETRA_CHK_ERR(-102);

//   double DN = N_;
//   double DNRHS = NRHS_;
//   if (!Solved()) EPETRA_CHK_ERR(-100); // Must have an existing solution
//   if (A_==AF_) EPETRA_CHK_ERR(-101); // Cannot apply refine if no original copy of A.

//   if (FERR_ != 0) delete [] FERR_; // Always start with a fresh copy of FERR_, since NRHS_ may change
//   FERR_ = new double[NRHS_];
//   if (BERR_ != 0) delete [] BERR_; // Always start with a fresh copy of BERR_, since NRHS_ may change
//   BERR_ = new double[NRHS_];
//   AllocateWORK();
//   AllocateIWORK();

//   LDB_ = LDX_ = N_;

//   std::cout << " ldb ldx "<<LDB_<<" "<<LDX_<<std::endl;

//   lapack.GERFS(TRANS_, N_, NRHS_, A_, LDA_, AF_, LDAF_, IPIV_,
// 	       B_, LDB_, X_, LDX_, FERR_, BERR_,
// 	       WORK_, IWORK_, &INFO_);

//   SolutionErrorsEstimated_ = true;
//   ReciprocalConditionEstimated_ = true;
//   SolutionRefined_ = true;

//   UpdateFlops(2.0*DN*DN*DNRHS); // Not sure of count

//   EPETRA_CHK_ERR(INFO_);
//   return(0);

 }

//=============================================================================
int Ifpack_SerialTriDiSolver::Invert(void)
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

  lapack.GETRI ( N_, AF_, LDAF_, IPIV_, WORK_, LWORK_, &INFO_);

  double DN = N_;
  UpdateFlops((DN*DN*DN));
  Inverted_ = true;
  Factored_ = false;

  EPETRA_CHK_ERR(INFO_);
  return(0);
}

//=============================================================================
int Ifpack_SerialTriDiSolver::ReciprocalConditionEstimate(double & Value)
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
  // We will assume a one-norm condition number \\ works for TriDi
  lapack.GECON( '1', N_, AF_, LDAF_, ANORM_, &RCOND_, WORK_, IWORK_, &INFO_);
  ReciprocalConditionEstimated_ = true;
  Value = RCOND_;
  UpdateFlops(2*N_*N_); // Not sure of count
  EPETRA_CHK_ERR(INFO_);
  return(0);
}
//=============================================================================
void Ifpack_SerialTriDiSolver::Print(std::ostream& os) const {

  if (Matrix_!=0) os << "Solver Matrix"          << std::endl << *Matrix_ << std::endl;
  if (Factor_!=0) os << "Solver Factored Matrix" << std::endl << *Factor_ << std::endl;
  if (LHS_   !=0) os << "Solver LHS"             << std::endl << *LHS_    << std::endl;
  if (RHS_   !=0) os << "Solver RHS"             << std::endl << *RHS_    << std::endl;

}
