
#include "Petra_BLAS_Easy_DLU_Problem.h"
//=============================================================================
Petra_BLAS_Easy_DLU_Problem::Petra_BLAS_Easy_DLU_Problem(void)
  : Petra_Flops(),
    Petra_BLAS(),
    Petra_LAPACK(),
    Transpose_(false),
    Factored_(false),
    Solved_(false),
    Inverted_(false),
    ReciprocalConditionEstimated_(false),
    
    LWORK_(0),
    
    IPIV_(0),
    IWORK_(0),
    WORK_(0),
    RCOND_(0.0),
    
    Operator_(0),
    RHS_(0),
    LHS_(0)

{
}

//=============================================================================
Petra_BLAS_Easy_DLU_Problem::Petra_BLAS_Easy_DLU_Problem(const Petra_BLAS_Easy_DLU_Problem& Source)
  : Petra_Flops(),  
    Petra_BLAS(),
    Petra_LAPACK(),
    Transpose_(Source.Transpose_),
    Factored_(Source.Factored_),
    Solved_(Source.Solved_),
    Inverted_(Source.Inverted_),
    ReciprocalConditionEstimated_(Source.ReciprocalConditionEstimated_),
    
    LWORK_(0),   
    IPIV_(0),
    IWORK_(0),
    WORK_(0),
    RCOND_(Source.RCOND_),
    
    Operator_(Source.Operator_),
    RHS_(Source.RHS_),
    LHS_(Source.LHS_)

{
}
//=============================================================================
Petra_BLAS_Easy_DLU_Problem::~Petra_BLAS_Easy_DLU_Problem()
{
  DeleteArrays();
}
//=============================================================================
void Petra_BLAS_Easy_DLU_Problem::DeleteArrays(void)
{
  if (IWORK_ != 0) {delete [] IWORK_; IWORK_ = 0;}
  if (IPIV_ != 0)  {delete [] IPIV_;IPIV_ = 0;}
  if (WORK_ != 0)  {delete [] WORK_;WORK_ = 0;}
}
//=============================================================================
int Petra_BLAS_Easy_DLU_Problem::SetOperator(Petra_BLAS_DGE_Matrix & A) {
  
  Operator_ = &A;
  Factored_ = false;
  Solved_ = false;
  return(0);
}

//=============================================================================
int Petra_BLAS_Easy_DLU_Problem::SetRHS(const Petra_BLAS_DGE_Matrix & B) {
  
  RHS_ = & (Petra_BLAS_DGE_Matrix &) B;
  Solved_ = false;
  return(0);
}

//=============================================================================
int Petra_BLAS_Easy_DLU_Problem::SetLHS(Petra_BLAS_DGE_Matrix & X) {
  
  LHS_ = &X;
  Solved_ = false;
  return(0);
}
//=============================================================================
int Petra_BLAS_Easy_DLU_Problem::Factor(void) {
  if (Factored()) return(0); // Already factored
  if (Operator_==0) PETRA_CHK_ERR(-1); // Matrix is not set
  if (Inverted()) PETRA_CHK_ERR(-100); // Cannot factor inverted matrix

  if (IPIV_==0) IPIV_ = new int[PETRA_MIN(Operator().M(),Operator().N())]; // Allocated Pivot vector if not already done.

  if (IPIV_==0) PETRA_CHK_ERR(-12);

  int M = Operator().M();
  int N = Operator().N();
  double * A = Operator().A();
  int LDA = Operator().LDA();
  int INFO = 0;
  GETRF (M, N, A, LDA, IPIV_, &INFO);

  Factored_ = true;
  double DM = M;
  double DN = N;
  UpdateFlops(DM*DN*DN - (DN*DN*DN)/3.0 - DN*DN/2.0 + 5.0*DN/6.0);

  PETRA_CHK_ERR(INFO);

}
//=============================================================================
int Petra_BLAS_Easy_DLU_Problem::Solve(void) {
  int ierr = 0;

  // We will call one of two routines depending on what services the user wants and 
  // whether or not the matrix has been inverted or factored already.
  //
  // If the matrix has been inverted, use DGEMM to compute solution.
  // Otherwise, if the matrix is already factored we will call the TRS interface.

  if (Operator_==0) PETRA_CHK_ERR(-2); // No A
  if (RHS_==0) PETRA_CHK_ERR(-3); // No B
  if (LHS_==0) PETRA_CHK_ERR(-4); // No X
  if (RHS().M()!=Operator().N() || RHS().N() != LHS().N()) PETRA_CHK_ERR(-5);
  if (Operator().A()==0) PETRA_CHK_ERR(-6);
  if (Operator().LDA()<1) PETRA_CHK_ERR(-7);
  if (RHS().A()==0) PETRA_CHK_ERR(-8);
  if (RHS().LDA()<1) PETRA_CHK_ERR(-9);
  if (LHS().A()==0) PETRA_CHK_ERR(-10);
  if (LHS().LDA()<1) PETRA_CHK_ERR(-11);

  int N = Operator().N();
  int NRHS = RHS().N();
  double DN = N;
  double DNRHS = NRHS;

  char TRANS = 'N'; if (Transpose()) TRANS = 'T';

  if (Inverted()) {

    if (RHS_==LHS_) PETRA_CHK_ERR(-100); // B and X must be different for this case

    GEMM(TRANS, 'N', N, NRHS, N, 1.0, Operator().A(), Operator().LDA(),
		RHS().A(), RHS().LDA(), 0.0, LHS().A(), LHS().LDA());
    UpdateFlops(2.0*DN*DN*DNRHS);
    Solved_ = true;
  }
  else {

    if (!Factored()) Factor(); // Matrix must be factored
    if (RHS_!=LHS_) *LHS_ = *RHS_;  
    
    int INFO = 0;
    GETRS(TRANS, N, NRHS, Operator().A(), Operator().LDA(), IPIV_, LHS().A(), LHS().LDA(), &INFO);
    if (INFO!=0) PETRA_CHK_ERR(INFO);
    UpdateFlops((2.0*DN*DN-DN)*DNRHS);
    Solved_ = true;

  }
  return(0);
}
//=============================================================================
int Petra_BLAS_Easy_DLU_Problem::Invert(void)
{
  if (!Factored()) Factor(); // Need matrix factored.

  /* This section works with LAPACK Version 3.0 only 
  // Setting LWORK = -1 and calling GETRI will return optimal work space size in WORK_TMP
  int LWORK_TMP = -1;
  double WORK_TMP;
  int INFO = 0;
  GETRI ( N_, AF_, LDAF_, IPIV_, &WORK_TMP, &LWORK_TMP, &INFO);
  LWORK_TMP = WORK_TMP; // Convert to integer
  if (LWORK_TMP>LWORK_) {
  if (WORK_!=0) delete WORK_;
  LWORK_ = LWORK_TMP;
  WORK_ = new double[LWORK_];
  }
  */
  // This section will work with any version of LAPACK 
  AllocateWORK();
  int INFO = 0;
  GETRI ( Operator().N(), Operator().A(), Operator().LDA(), IPIV_, WORK_, &LWORK_, &INFO);

  double DN = Operator().N();
  UpdateFlops((DN*DN*DN));
  Inverted_ = true;
  Factored_ = false;
  
  PETRA_CHK_ERR(INFO);
}

//=============================================================================
int Petra_BLAS_Easy_DLU_Problem::ReciprocalConditionEstimate(double & Value)
{
  int ierr = 0;
  if (ReciprocalConditionEstimated()) {
    Value = RCOND_;
    return(0); // Already computed, just return it.
  }

  if (!Factored()) ierr = Factor(); // Need matrix factored.
  if (ierr!=0) PETRA_CHK_ERR(ierr-2);
  double anorm = Operator().OneNorm();

  AllocateWORK();
  AllocateIWORK();
  // We will assume a one-norm condition number
  int INFO = 0;
  GECON( '1', Operator().N(), Operator().A(), Operator().LDA(), anorm, &RCOND_, WORK_, IWORK_, &INFO);
  ReciprocalConditionEstimated_ = true;
  Value = RCOND_;
  int N = Operator().N();
  UpdateFlops(2*N*N); // Not sure of count
  PETRA_CHK_ERR(INFO);
}

//=============================================================================
void Petra_BLAS_Easy_DLU_Problem::Print (ostream& os) const
{
  os << "Linear Operator  A: " << *Operator_ << endl << endl;
  os << "Left-hand-side   X: " << *LHS_ << endl << endl;
  os << "Right-hand-side  B: " << *RHS_ << endl << endl;
  return;
}
