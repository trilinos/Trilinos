
#include "Petra_BLAS_DGE_Matrix.h"
//=============================================================================
Petra_BLAS_DGE_Matrix::Petra_BLAS_DGE_Matrix(void)
  : Petra_Flops(),
    M_(0),
    N_(0),
    LDA_(0),
    A_Copied_(false),
    A_(0)
{
}

//=============================================================================
Petra_BLAS_DGE_Matrix::Petra_BLAS_DGE_Matrix(Petra_DataAccess CV, double *A, int LDA, 
					     int NumRows, int NumCols)
  : Petra_Flops(),
    M_(NumRows),
    N_(NumCols),
    LDA_(LDA),
    A_Copied_(false),    
    A_(A)

{
  if (CV==Copy) {
    LDA_ = M_;
    A_ = new double[LDA_*N_];
    CopyMat(A, LDA, M_, N_, A_, LDA_);
    A_Copied_ = true;
  }

}
//=============================================================================
Petra_BLAS_DGE_Matrix::Petra_BLAS_DGE_Matrix(const Petra_BLAS_DGE_Matrix& Source)
  : Petra_Flops(),  
    M_(Source.M_),
    N_(Source.N_),
    LDA_(Source.LDA_),
    A_Copied_(true),
    A_(Source.A_)

{

  int i;
  LDA_ = M_;
  A_ = new double[LDA_*N_];
  CopyMat(Source.A_, Source.LDA_, M_, N_, A_, LDA_);
}
int Petra_BLAS_DGE_Matrix::Reshape(int NumRows, int NumCols) {

  // Allocate space for new matrix
  double * A_tmp = new double[NumRows*NumCols];
  for (int k = 0; k < NumRows*NumCols; k++) A_tmp[k] = 0.0; // Zero out values
  int M_tmp = PETRA_MIN(M_, NumRows);
  int N_tmp = PETRA_MIN(N_, NumCols);
  if (A_ != 0) CopyMat(A_, LDA_, M_tmp, N_tmp, A_tmp, NumRows); // Copy principal submatrix of A to new A
  
  DeleteArrays(); // Get rid of anything that might be already allocated
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
  A_ = A_tmp; // Set pointer to new A
  A_Copied_ = true;

  return(0);
}
int Petra_BLAS_DGE_Matrix::Shape(int NumRows, int NumCols) {
  DeleteArrays(); // Get rid of anything that might be already allocated
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
  A_ = new double[LDA_*N_];
  for (int k = 0; k < LDA_*N_; k++) A_[k] = 0.0; // Zero out values
  A_Copied_ = true;

  return(0);
}
//=============================================================================
Petra_BLAS_DGE_Matrix::~Petra_BLAS_DGE_Matrix()
{
  DeleteArrays();
}
//=============================================================================
void Petra_BLAS_DGE_Matrix::DeleteArrays(void)
{
  if (A_Copied_)   {delete [] A_; A_ = 0; A_Copied_ = false;}
}
//=============================================================================
void Petra_BLAS_DGE_Matrix::CopyMat(double * A, int LDA, int NumRows, int NumCols, 
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
double Petra_BLAS_DGE_Matrix::OneNorm(void) {

  int i, j;

    double anorm = 0.0;
    double * ptr;
    for (j=0; j<N_; j++) {
      double sum=0.0;
      ptr = A_ + j*LDA_;
      for (i=0; i<M_; i++) sum += fabs(*ptr++);
      anorm = PETRA_MAX(anorm, sum);
    }
    UpdateFlops(N_*N_);
    return(anorm);
}
//=============================================================================
double Petra_BLAS_DGE_Matrix::InfNorm(void) {

  int i, j;

    double anorm = 0.0;
    double * ptr;

    // Loop across columns in inner loop.  Most expensive memory access, but 
    // requires no extra storage.
    for (i=0; i<M_; i++) {
      double sum=0.0;
      ptr = A_ + i;
      for (j=0; j<N_; j++) {
	sum += fabs(*ptr);
	ptr += LDA_;
      }
      anorm = PETRA_MAX(anorm, sum);
    }
    UpdateFlops(N_*N_);
    return(anorm);
}
//=========================================================================
double& Petra_BLAS_DGE_Matrix::operator () (int RowIndex, int ColIndex)  {

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
const double& Petra_BLAS_DGE_Matrix::operator () (int RowIndex, int ColIndex) const  {

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
const double* Petra_BLAS_DGE_Matrix::operator [] (int ColIndex) const  {

  if (ColIndex>=N_) {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << N_-1;
    abort();
  }
   return(A_ + ColIndex*LDA_);
}
//=========================================================================
double* Petra_BLAS_DGE_Matrix::operator [] (int ColIndex)  {

  if (ColIndex>=N_) {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << N_-1;
    abort();
  }
   return(A_+ ColIndex*LDA_);
}
//=========================================================================
int  Petra_BLAS_DGE_Matrix::Multiply (char TransA, char TransB, double ScalarAB, 
				      const Petra_BLAS_DGE_Matrix& A, 
				      const Petra_BLAS_DGE_Matrix& B,
				      double Scalar ) {
  // Check for compatible dimensions
  
  int A_nrows = (TransA=='T') ? A.N() : A.M();
  int A_ncols = (TransA=='T') ? A.M() : A.N();
  int B_nrows = (TransB=='T') ? B.N() : B.M();
  int B_ncols = (TransB=='T') ? B.M() : B.N();
  
  if (M_        != A_nrows     ||
      A_ncols   != B_nrows     ||
      N_        != B_ncols  ) PETRA_CHK_ERR(-1); // Return error

    
  // Call GEMM function
  GEMM(TransA, TransB, M_, N_, A_ncols, ScalarAB, A.A(), A.LDA(), 
       B.A(), B.LDA(), Scalar, A_, LDA_);
  long int nflops = 2*M_;
  nflops *= N_;
  nflops *= A_ncols;
  UpdateFlops(nflops);

  return(0);
}
void Petra_BLAS_DGE_Matrix::Print(ostream& os) const {

  for (int i=0; i<M_; i++) {
    for (int j=0; j<N_; j++){
      os << (*this)(i,j);
    }
    os << endl;
  }
  return;
}
