
#include "Petra_RDP_SPD_DenseMatrix.h"
//=============================================================================
Petra_RDP_SPD_DenseMatrix::Petra_RDP_SPD_DenseMatrix(void)
  : Petra_RDP_DenseMatrix(),
    Upper_(false),
    UPLO_('L'),
    SCOND_(-1.0)

{
}
//=============================================================================
Petra_RDP_SPD_DenseMatrix::Petra_RDP_SPD_DenseMatrix(Petra_DataAccess CV, double *A, int LDA, int NumRowsCols)
  : Petra_RDP_DenseMatrix(CV, A, LDA, NumRowsCols, NumRowsCols),
    Upper_(false),
    UPLO_('L'),
    SCOND_(-1.0)

{
}
//=============================================================================
Petra_RDP_SPD_DenseMatrix::Petra_RDP_SPD_DenseMatrix(const Petra_RDP_SPD_DenseMatrix& Source)
  : Petra_RDP_DenseMatrix(Source),  
    Upper_(Source.Upper_),
    UPLO_(Source.UPLO_),
    SCOND_(Source.SCOND_)
{
}
//=============================================================================
Petra_RDP_SPD_DenseMatrix::~Petra_RDP_SPD_DenseMatrix()
{
}
//=============================================================================
int Petra_RDP_SPD_DenseMatrix::Factor(void) {
  if (Factored()) return(0); // Already factored
  if (Inverted()) return(-100); // Cannot factor inverted matrix
  int ierr = 0;

  ierr = OneNorm();
  if (ierr!=0) return(ierr-1); // Compute 1-Norm of A


  // If we want to refine the solution, then the factor must
  // be stored separatedly from the original matrix

  if (A_ == AF_)
    if (RefineSolution_ ) {
      LDAF_ = N_;
      AF_ = new double[LDAF_*N_];
      CopyMat(A_, LDA_, N_, N_, AF_, LDAF_);
    }
  if (Equilibrate_) ierr = Equilibrate_A();

  if (ierr!=0) return(ierr-2);
  
  POTRF (UPLO_, N_, AF_, LDAF_, &INFO_);
  Factored_ = true;
  double DN = N_;
  UpdateFlops((DN*DN*DN)/3.0);

  return(INFO_);

}

//=============================================================================
int Petra_RDP_SPD_DenseMatrix::Solve(void) {
  int ierr = 0;

  // We will call one of four routines depending on what services the user wants and 
  // whether or not the matrix has been inverted or factored already.
  //
  // If the matrix has been inverted, use DGEMM to compute solution.
  // Otherwise, if the user want the matrix to be equilibrated or wants a refined solution, we will
  // call the X interface.
  // Otherwise, if the matrix is already factored we will call the TRS interface.
  // Otherwise, if the matrix is unfactored we will call the SV interface.

  double DN = N_;
  double DNRHS = NRHS_;
  if (A_Equilibrated_ && !B_Equilibrated_) return(-1); // Matrix and vectors must be similarly scaled
  if (!A_Equilibrated_ && B_Equilibrated_) return(-2);
  if (B_==0) return(-3); // No B
  if (X_==0) return(-4); // No B

  if (ShouldEquilibrate() && !A_Equilibrated_) ierr = 1; // Warn that the system should be equilibrated.

  if (Inverted()) {

    if (B_==X_) return(-100); // B and X must be different for this case
    GEMM('N', 'N', N_, NRHS_, N_, 1.0, AF_, LDAF_,
		B_, LDB_, 0.0, X_, LDX_);
    if (INFO_!=0) return(INFO_);
    UpdateFlops(2.0*DN*DN*DNRHS);
    Solved_ = true;
  }
  else {

    if (!Factored()) Factor(); // Matrix must be factored
    if (B_!=X_) CopyMat(B_, LDB_, N_, NRHS_, X_, LDX_);  
    
    POTRS(UPLO_, N_, NRHS_, AF_, LDAF_, X_, LDX_, &INFO_);
    if (INFO_!=0) return(INFO_);
    UpdateFlops(2.0*DN*DN*DNRHS);
    Solved_ = true;

  }
  int ierr1=0;
  if (RefineSolution_) ierr1 = ApplyRefinement();
  if (ierr1!=0) return(ierr1);
  else
    return(ierr);
}
//=============================================================================
int Petra_RDP_SPD_DenseMatrix::ApplyRefinement(void)
{
  double DN = N_;
  double DNRHS = NRHS_;
  if (!Solved()) return(-100); // Must have an existing solution
  if (A_==AF_) return(-101); // Cannot apply refine if no original copy of A.

  if (FERR_ != 0) delete FERR_; // Always start with a fresh copy of FERR_, since NRHS_ may change
  FERR_ = new double[NRHS_];
  if (BERR_ != 0) delete BERR_; // Always start with a fresh copy of FERR_, since NRHS_ may change
  BERR_ = new double[NRHS_];
  AllocateWORK();
  AllocateIWORK();
  
  PORFS(UPLO_, N_, NRHS_, A_, LDA_, AF_, LDAF_,
	       B_, LDB_, X_, LDX_, FERR_, BERR_, 
	       WORK_, IWORK_, &INFO_);
  
  
  SolutionErrorsEstimated_ = true;
  ReciprocalConditionEstimated_ = true;
  SolutionRefined_ = true;
  
  UpdateFlops(2.0*DN*DN*DNRHS); // Not sure of count
  
  return(INFO_);
  
}

//=============================================================================
int Petra_RDP_SPD_DenseMatrix::ComputeEquilibrateScaling(void) {
  if (R_!=0) return(0); // Already computed
 
  double DN = N_;
  R_ = new double[N_];
  C_ = R_;
  
  POEQU (N_, AF_, LDAF_, R_, &SCOND_, &AMAX_, &INFO_);
  if (INFO_ != 0) return(INFO_);

  if (SCOND_<0.1 || AMAX_ < Underflow_ || AMAX_ > Overflow_) ShouldEquilibrate_ = true;

  UpdateFlops(2.0*DN*DN);
  
  return(0);
}

//=============================================================================
int Petra_RDP_SPD_DenseMatrix::Equilibrate_A(void)
{
  int i, j;
  int ierr = 0;

  double DN = N_;

  if (A_Equilibrated_) return(0); // Already done
  if (R_==0) ierr = ComputeEquilibrateScaling(); // Compute S if needed
  if (ierr!=0) return(ierr);
  if (A_==AF_) {
    double * ptr;
    for (j=0; j<N_; j++) {
      ptr = A_ + j*LDA_;
      double s1 = R_[j];
      for (i=0; i<N_; i++) {
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
      for (i=0; i<N_; i++) {
	*ptr = *ptr*s1*R_[i];
	ptr++;
	*ptr1 = *ptr1*s1*R_[i];
	ptr1++;
      }
    }
  }
  
  A_Equilibrated_ = true;
  UpdateFlops(2.0*DN*DN);
  
  return(0);
}

//=============================================================================
int Petra_RDP_SPD_DenseMatrix::Equilibrate_B(void)
{
  int i, j;
  int ierr = 0;

  if (B_Equilibrated_) return(0); // Already done
  if (R_==0) ierr = ComputeEquilibrateScaling(); // Compute S if needed
  if (ierr!=0) return(ierr);

  double * ptr;
  for (j=0; j<NRHS_; j++) {
    ptr = B_ + j*LDB_;
    for (i=0; i<N_; i++) {
      *ptr = *ptr*R_[i];
      ptr++;
    }
  }

  
  B_Equilibrated_ = true;
  UpdateFlops((double) N_ * (double) NRHS_);
  
  return(0);
}

//=============================================================================
int Petra_RDP_SPD_DenseMatrix::Unequilibrate_X(void)
{
  int i, j;

  if (!B_Equilibrated_) return(0); // Nothing to do

  double * ptr;
  for (j=0; j<NRHS_; j++) {
    ptr = X_ + j*LDX_;
    for (i=0; i<N_; i++) {
      *ptr = *ptr*R_[i];
      ptr++;
    }
  }

  
  UpdateFlops((double) N_ * (double) NRHS_);
  
  return(0);
}

//=============================================================================
int Petra_RDP_SPD_DenseMatrix::Invert(void)
{
  if (!Factored()) Factor(); // Need matrix factored.
  POTRI ( UPLO_, N_, AF_, LDAF_, &INFO_);
  CopyUPLOMat(Upper_, AF_, LDAF_, N_);
  double DN = N_;
  UpdateFlops((DN*DN*DN));
  Inverted_ = true;
  Factored_ = false;
  
  return(INFO_);
}

//=============================================================================
int Petra_RDP_SPD_DenseMatrix::ReciprocalConditionEstimate(double & Value)
{
  int ierr = 0;
  if (ReciprocalConditionEstimated()) {
    Value = RCOND_;
    return(0); // Already computed, just return it.
  }

  if (ANORM_<0.0) ierr = OneNorm();
  if (ierr!=0) return(ierr-1);
  if (!Factored()) ierr = Factor(); // Need matrix factored.
  if (ierr!=0) return(ierr-2);

  AllocateWORK();
  AllocateIWORK();
  POCON( UPLO_, N_, AF_, LDAF_, &ANORM_, &RCOND_, WORK_, IWORK_, &INFO_);
  ReciprocalConditionEstimated_ = true;
  Value = RCOND_;
  UpdateFlops(2*N_*N_); // Not sure of count
  return(INFO_);
}
//=============================================================================
void Petra_RDP_SPD_DenseMatrix::CopyUPLOMat(bool Upper, double * A, int LDA, int NumRows) {

  int i, j;
  double * ptr1;
  double * ptr2;

  if (Upper) {
    for (j=1; j<NumRows; j++) {
      ptr1 = A + j;
      ptr2 = A + j*LDA;
      for (i=0; i<j; i++) {
	*ptr1 = *ptr2++;
	ptr1+=LDA;
      }
    }
  }
  else {
    for (i=1; i<NumRows; i++) {
      ptr1 = A + i;
      ptr2 = A + i*LDA;
      for (j=0; j<i; j++) {
	*ptr2++ = *ptr1;
	ptr1+=LDA;
      }
    }
  }
}

