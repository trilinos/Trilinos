
#include "Petra_RDP_DenseMatrix.h"
//=============================================================================
Petra_RDP_DenseMatrix::Petra_RDP_DenseMatrix(void)
  : Petra_Flops(),
    Petra_BLAS(),
    Petra_LAPACK(),
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
    A_Copied_(false),
    B_Copied_(false),
    X_Copied_(false),
    
    TRANS_('N'),
    
    M_(0),
    N_(0),
    Min_MN_(minfn(M_,N_)),
    NRHS_(0),
    LDA_(0),
    LDAF_(0),
    LDB_(0),
    LDX_(0),
    INFO_(0),
    LWORK_(0),
    
    IPIV_(0),
    IWORK_(0),
    
    ANORM_(-1.0),
    RCOND_(-1.0),
    ROWCND_(-1.0),
    COLCND_(-1.0),
    AMAX_(-1.0),
    
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
}

//=============================================================================
Petra_RDP_DenseMatrix::Petra_RDP_DenseMatrix(Petra_DataAccess CV, double *A, int LDA, 
					     int NumRows, int NumCols)
  : Petra_Flops(),
    Petra_BLAS(),
    Petra_LAPACK(),
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
    A_Copied_(false),
    B_Copied_(false),
    X_Copied_(false),
    
    TRANS_('N'),
    
    M_(NumRows),
    N_(NumCols),
    Min_MN_(minfn(M_,N_)),
    NRHS_(0),
    LDA_(LDA),
    LDAF_(LDA),
    LDB_(0),
    LDX_(0),
    INFO_(0),
    LWORK_(0),
    
    IPIV_(0),
    IWORK_(0),
    
    ANORM_(-1.0),
    RCOND_(-1.0),
    ROWCND_(-1.0),
    COLCND_(-1.0),
    AMAX_(-1.0),
    
    A_(A),
    FERR_(0),
    BERR_(0),
    AF_(A),
    WORK_(0),
    R_(0),
    C_(0),
    B_(0),
    X_(0)

{
  if (CV==Copy) {
    LDA_ = M_;
    A_ = new double[LDA_*N_];
    CopyMat(A, LDA, M_, N_, A_, LDA_);
    AF_ = A_;
    LDAF_ = LDA_;
    A_Copied_ = true;
  }

}
//=============================================================================
Petra_RDP_DenseMatrix::Petra_RDP_DenseMatrix(const Petra_RDP_DenseMatrix& Source)
  : Petra_Flops(),  
    Petra_BLAS(),
    Petra_LAPACK(),
    Equilibrate_(Source.Equilibrate_),
    ShouldEquilibrate_(Source.ShouldEquilibrate_),
    A_Equilibrated_(Source.A_Equilibrated_),
    B_Equilibrated_(Source.B_Equilibrated_),
    Transpose_(Source.Transpose_),
    Factored_(Source.Factored_),
    EstimateSolutionErrors_(Source.EstimateSolutionErrors_),
    SolutionErrorsEstimated_(Source.SolutionErrorsEstimated_),
    Solved_(Source.Solved_),
    Inverted_(Source.Inverted_),
    ReciprocalConditionEstimated_(Source.ReciprocalConditionEstimated_),
    RefineSolution_(Source.RefineSolution_),
    SolutionRefined_(Source.SolutionRefined_),
    A_Copied_(true),
    B_Copied_(Source.B_Copied_),
    X_Copied_(Source.X_Copied_),
    
    TRANS_(Source.TRANS_),

    M_(Source.M_),
    N_(Source.N_),
    Min_MN_(Source.Min_MN_),
    NRHS_(0),
    LDA_(Source.LDA_),
    LDAF_(Source.LDAF_),
    LDB_(0),
    LDX_(0),
    INFO_(0),
    LWORK_(0),
    
    IPIV_(0),
    IWORK_(0),
    
    ANORM_(Source.ANORM_),
    RCOND_(Source.RCOND_),
    ROWCND_(Source.ROWCND_),
    COLCND_(Source.COLCND_),
    AMAX_(Source.AMAX_),
    
    A_(Source.A_),
    FERR_(0),
    BERR_(0),
    AF_(Source.AF_),
    WORK_(0),
    R_(Source.R_),
    C_(Source.C_),
    B_(Source.B_),
    X_(Source.X_)

{

  int i;
  LDA_ = M_;
  A_ = new double[LDA_*N_];
  CopyMat(Source.A_, Source.LDA_, M_, N_, A_, LDA_);
  if (Source.AF_ !=Source.A_) {
    LDAF_ = M_;
    AF_ = new double[LDAF_*N_];
    CopyMat(Source.AF_, Source.LDAF_, M_, N_, AF_, LDAF_);
  }
  else { AF_ = A_; LDAF_ = LDA_;}

  if (IPIV_ != 0) {
    IPIV_ = new int[Min_MN_];
    for (i=0; i<Min_MN_; i++) IPIV_[i] = Source.IPIV_[i];
  }
  if (R_ != 0) {
    R_ = new double[M_];
    CopyMat(Source.R_, M_, M_, 1, R_, M_);
  }
  if (C_ != 0) {
    C_ = new double[N_];
    CopyMat(Source.C_, N_, N_, 1, C_, N_);
  }
  if (B_ != 0) {
    B_ = new double[N_*NRHS_];
    LDB_ = N_;
    CopyMat(Source.B_, Source.LDB_, N_, NRHS_, B_, LDB_);
    B_Copied_ = true;
  }
  if (X_ != 0) {
    if (Source.B_==Source.X_) {  // If source has X and B shared, copy will also
      X_ = B_;
      LDX_ = LDB_;
    }
    else {
      X_ = new double[M_*NRHS_];
      LDX_ = M_;
      CopyMat(Source.X_, Source.LDX_, M_, NRHS_, X_, LDX_);
      X_Copied_ = true;
    }
  }
}
int Petra_RDP_DenseMatrix::Reshape(int NumRows, int NumCols) {

  // Allocate space for new matrix
  double * A_tmp = new double[NumRows*NumCols];
  for (int k = 0; k < NumRows*NumCols; k++) A_tmp[k] = 0.0; // Zero out values
  int M_tmp = minfn(M_, NumRows);
  int N_tmp = minfn(N_, NumCols);
  if (A_ != 0) CopyMat(A_, LDA_, M_tmp, N_tmp, A_tmp, NumRows); // Copy principal submatrix of A to new A
  
  DeleteArrays(); // Get rid of anything that might be already allocated
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
  Min_MN_ = minfn(M_, N_);
  A_ = A_tmp; // Set pointer to new A
  AF_ = A_;
  LDAF_ = LDA_;
  A_Copied_ = true;

  return(0);
}
int Petra_RDP_DenseMatrix::Shape(int NumRows, int NumCols) {
  DeleteArrays(); // Get rid of anything that might be already allocated
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
  Min_MN_ = minfn(M_, N_);
  A_ = new double[LDA_*N_];
  for (int k = 0; k < LDA_*N_; k++) A_[k] = 0.0; // Zero out values
  AF_ = A_;
  LDAF_ = LDA_;
  A_Copied_ = true;

  return(0);
}
//=============================================================================
Petra_RDP_DenseMatrix::~Petra_RDP_DenseMatrix()
{
  DeleteArrays();
}
//=============================================================================
void Petra_RDP_DenseMatrix::DeleteArrays(void)
{
  if (IWORK_ != 0) {delete [] IWORK_; IWORK_ = 0;}
  if (FERR_ != 0)  {delete [] FERR_; FERR_ = 0;}
  if (BERR_ != 0)  {delete [] BERR_; BERR_ = 0;}
  if (AF_ != A_ && AF_ != 0)   {delete [] AF_; AF_ = 0;}
  if (AF_ !=0) AF_ = 0;
  if (A_Copied_)   {delete [] A_; A_ = 0; A_Copied_ = false;}
  if (IPIV_ != 0)  {delete [] IPIV_;IPIV_ = 0;}
  if (WORK_ != 0)  {delete [] WORK_;WORK_ = 0;}
  if (R_ != 0 && R_ != C_)     {delete [] R_;R_ = 0;}
  if (R_ != 0) R_ = 0;
  if (C_ != 0)     {delete [] C_;C_ = 0;}
  if (B_Copied_)   {delete [] B_; B_ = 0; B_Copied_ = false;}
  if (X_Copied_)   {delete [] X_; X_ = 0; X_Copied_ = false;}
}
//=============================================================================
int Petra_RDP_DenseMatrix::Factor(void) {
  if (Factored()) return(0); // Already factored
  if (Inverted()) return(-100); // Cannot factor inverted matrix
  int ierr = 0;

  ierr = OneNorm();
  if (ierr!=0) return(ierr-1); // Compute 1-Norm of A


  // If we want to refine the solution, then the factor must
  // be stored separatedly from the original matrix

  if (A_ == AF_)
    if (RefineSolution_ ) {
      LDAF_ = M_;
      AF_ = new double[LDAF_*N_];
      CopyMat(A_, LDA_, M_, N_, AF_, LDAF_);
    }
  
  if (Equilibrate_) ierr = Equilibrate_A();

  if (ierr!=0) return(ierr-2);
  
  if (IPIV_==0) IPIV_ = new int[Min_MN_]; // Allocated Pivot vector if not already done.

  GETRF (M_, N_, AF_, LDAF_, IPIV_, &INFO_);

  Factored_ = true;
  double DN = N_;
  UpdateFlops(2.0*(DN*DN*DN)/3.0);

  return(INFO_);

}

//=============================================================================
int Petra_RDP_DenseMatrix::SetVectors(const Petra_RDP_DenseMatrix & B, Petra_RDP_DenseMatrix & X) {
  int ierr = 0;

  
  if (B.M()!=X.M() || B.N() != X.N()) return(-1);
  if (B.A()==0) return(-2);
  if (B.LDA()<1) return(-3);
  if (X.A()==0) return(-4);
  if (X.LDA()<1) return(-5);

  if (B_Copied_)   delete [] B_; // delete old B if any
  if (X_Copied_)   delete [] X_; // delete old X if any

  NRHS_ = B.N();

  B_ = B.A();
  LDB_ = B.LDA();
  B_Copied_ = false;
  X_ = X.A();
  LDX_ = X.LDA();
  X_Copied_ = false;

  Solved_ = false;
  SolutionErrorsEstimated_ = false;
  SolutionRefined_ = false;
  B_Equilibrated_ = false;

  if (Equilibrate_) {
    ierr = Equilibrate_B();
    B_Equilibrated_ = true;
  }
  return(ierr);
}

//=============================================================================
int Petra_RDP_DenseMatrix::Solve(void) {
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

    GEMM(TRANS_, 'N', N_, NRHS_, N_, 1.0, AF_, LDAF_,
		B_, LDB_, 0.0, X_, LDX_);
    if (INFO_!=0) return(INFO_);
    UpdateFlops(2.0*DN*DN*DNRHS);
    Solved_ = true;
  }
  else {

    if (!Factored()) Factor(); // Matrix must be factored
    if (B_!=X_) CopyMat(B_, LDB_, N_, NRHS_, X_, LDX_);  
    
    GETRS(TRANS_, N_, NRHS_, AF_, LDAF_, IPIV_, X_, LDX_, &INFO_);
    if (INFO_!=0) return(INFO_);
    UpdateFlops(2.0*DN*DN*DNRHS);
    Solved_ = true;

  }
  int ierr1=0;
  if (RefineSolution_ && !Inverted()) ierr1 = ApplyRefinement();
  if (ierr1!=0) return(ierr1);
  else
    return(ierr);
}
//=============================================================================
int Petra_RDP_DenseMatrix::ApplyRefinement(void)
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
  
  GERFS(TRANS_, N_, NRHS_, A_, LDA_, AF_, LDAF_, IPIV_,
	       B_, LDB_, X_, LDX_, FERR_, BERR_, 
	       WORK_, IWORK_, &INFO_);
  
  
  SolutionErrorsEstimated_ = true;
  ReciprocalConditionEstimated_ = true;
  SolutionRefined_ = true;
  
  UpdateFlops(2.0*DN*DN*DNRHS); // Not sure of count
  
  return(INFO_);
  
}

//=============================================================================
int Petra_RDP_DenseMatrix::ComputeEquilibrateScaling(void) {
  if (R_!=0) return(0); // Already computed
 
  double DM = M_;
  double DN = N_;
  R_ = new double[M_];
  C_ = new double[N_];
  
  GEEQU (M_, N_, AF_, LDAF_, R_, C_, &ROWCND_, &COLCND_, &AMAX_, &INFO_);
  if (INFO_ != 0) return(INFO_);

  if (COLCND_<0.1 || ROWCND_<0.1 || AMAX_ < Underflow_ || AMAX_ > Overflow_) ShouldEquilibrate_ = true;

  UpdateFlops(4.0*DM*DN);
  
  return(0);
}

//=============================================================================
int Petra_RDP_DenseMatrix::Equilibrate_A(void)
{
  int i, j;
  int ierr = 0;

  double DN = N_;
  double DM = M_;

  if (A_Equilibrated_) return(0); // Already done
  if (R_==0) ierr = ComputeEquilibrateScaling(); // Compute R and C if needed
  if (ierr!=0) return(ierr);
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
int Petra_RDP_DenseMatrix::Equilibrate_B(void)
{
  int i, j;
  int ierr = 0;

  if (B_Equilibrated_) return(0); // Already done
  if (R_==0) ierr = ComputeEquilibrateScaling(); // Compute R and C if needed
  if (ierr!=0) return(ierr);

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
int Petra_RDP_DenseMatrix::Unequilibrate_X(void)
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
int Petra_RDP_DenseMatrix::Invert(void)
{
  if (!Factored()) Factor(); // Need matrix factored.

  /* This section work with LAPACK Version 3.0 only 
  // Setting LWORK = -1 and calling GETRI will return optimal work space size in WORK_TMP
  int LWORK_TMP = -1;
  double WORK_TMP;
  GETRI ( N_, AF_, LDAF_, IPIV_, &WORK_TMP, &LWORK_TMP, &INFO_);
  LWORK_TMP = WORK_TMP; // Convert to integer
  if (LWORK_TMP>LWORK_) {
  if (WORK_!=0) delete WORK_;
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
  
  return(INFO_);
}

//=============================================================================
int Petra_RDP_DenseMatrix::ReciprocalConditionEstimate(double & Value)
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
  // We will assume a one-norm condition number
  GECON( '1', N_, AF_, LDAF_, ANORM_, &RCOND_, WORK_, IWORK_, &INFO_);
  ReciprocalConditionEstimated_ = true;
  Value = RCOND_;
  UpdateFlops(2*N_*N_); // Not sure of count
  return(INFO_);
}
//=============================================================================
void Petra_RDP_DenseMatrix::CopyMat(double * A, int LDA, int NumRows, int NumCols, 
					double * B, int LDB) {

  int i, j;
  double * ptr1 = B;
  double * ptr2;
  for (j=0; j<NumCols; j++) {
    ptr1 = B + j*LDB;
    ptr2 = A + j*LDA;
    for (i=0; i<NumRows; i++) *ptr1++ = *ptr2++;
  }
  return;
}
//=============================================================================
int Petra_RDP_DenseMatrix::OneNorm(void) {

  int i, j;

  if (!Factored() || A_ != AF_) { // A has unfactored values

    ANORM_ = 0.0;
    double * ptr;
    for (j=0; j<N_; j++) {
      double sum=0.0;
      ptr = A_ + j*LDA_;
      for (i=0; i<M_; i++) sum += fabs(*ptr++);
      ANORM_ = maxfn(ANORM_, sum);
    }
    UpdateFlops(N_*N_);
    return(0);
  }
  
  return(-100);
}
//=========================================================================
double& Petra_RDP_DenseMatrix::operator () (int RowIndex, int ColIndex)  {

  if (RowIndex>=M_) {
    cout << "Row index = " << RowIndex << " Out of Range 0 - " << M_-1;
    abort();
  }
  if (ColIndex>=N_) {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << N_-1;
    abort();
  }
   return(A_[ColIndex*LDA_ + RowIndex]);
}

//=========================================================================
const double& Petra_RDP_DenseMatrix::operator () (int RowIndex, int ColIndex) const  {

  if (RowIndex>=M_) {
    cout << "Row index = " << RowIndex << " Out of Range 0 - " << M_-1;
    abort();
  }
  if (ColIndex>=N_) {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << N_-1;
    abort();
  }
   return(A_[ColIndex*LDA_ + RowIndex]);
}
//=========================================================================
const double* Petra_RDP_DenseMatrix::operator [] (int ColIndex) const  {

  if (ColIndex>=N_) {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << N_-1;
    abort();
  }
   return(A_ + ColIndex*LDA_);
}
//=========================================================================
double* Petra_RDP_DenseMatrix::operator [] (int ColIndex)  {

  if (ColIndex>=N_) {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << N_-1;
    abort();
  }
   return(A_+ ColIndex*LDA_);
}
//=========================================================================
int  Petra_RDP_DenseMatrix::Multiply (char TransA, char TransB, double ScalarAB, 
				      const Petra_RDP_DenseMatrix& A, 
				      const Petra_RDP_DenseMatrix& B,
				      double Scalar ) {
  // Check for compatible dimensions
  
  int A_nrows = (TransA=='T') ? A.N() : A.M();
  int A_ncols = (TransA=='T') ? A.M() : A.N();
  int B_nrows = (TransB=='T') ? B.N() : B.M();
  int B_ncols = (TransB=='T') ? B.M() : B.N();
  
  if (M_        != A_nrows     ||
      A_ncols   != B_nrows     ||
      N_        != B_ncols  ) return(-1); // Return error
  
  // Call GEMM function
  GEMM(TransA, TransB, M_, N_, A_ncols, ScalarAB, A.A(), A.LDA(), 
       B.A(), B.LDA(), Scalar, A_, LDA_);
  double DM = M_;
  double DN = N_;
  double DK = A_ncols;
  UpdateFlops(2.0*DM*DN*DK);

  return(0);
}

// Non-member functions

ostream& operator << (ostream& os, const Petra_RDP_DenseMatrix& A)
{
  const int M = A.M();
  const int N = A.N();
  long olda = os.setf(ios::right,ios::adjustfield);
  long oldf = os.setf(ios::scientific,ios::floatfield);
  int oldp = os.precision(12);

  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++){
      os.width(20);
      os << A(i,j);
    }
    os << endl;
  }

  // Reset os flags

  os.setf(olda,ios::adjustfield);
  os.setf(oldf,ios::floatfield);
  os.precision(oldp);

  return os;
}
