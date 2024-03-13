#include "matrix.hpp"
#include <stk_util/util/BlasLapack.hpp>

#include <algorithm>
#include <cstring>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

void matvec(const double alpha, const Matrix<double>& a, const double* x, const double beta, double* y,
            BlasTrans transIn)
{
  char reverseTransBecauseAisRowMajor = transIn==BlasTrans::NoTrans ? 'T' : 'N';
  int m = a.extent0();
  int n = a.extent1();
  std::swap(m,n);
  int lda               = a.extent1();
  int incx              = 1;
  int incy              = 1;
  SIERRA_FORTRAN(dgemv)(&reverseTransBecauseAisRowMajor, &m, &n, &alpha, a.data(), &lda, x, &incx, &beta, y, &incy);
}

void matvec(BlasTrans trans, int m, int n, double alpha, const double* a, const double* x, double beta, double* y)
{
  char transChar = trans == BlasTrans::NoTrans ? 'N' : 'T';
  int lda = m;
  int incx = 1;
  int incy = 1;
  SIERRA_FORTRAN(dgemv)(&transChar, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

void solve_least_squares(Matrix<double>& a, Matrix<double>& b, Matrix<double>& work, LapackTrans trans)
{
  char transC = static_cast<char>(trans);
#ifndef NDEBUG
  int maxDim = std::max(a.extent0(), a.extent1());
  assert(b.extent0() >= maxDim);
  assert(b.extent1() == 1);
#endif
  Matrix<double> at = transpose(a);
  int m             = a.extent0();
  int n             = a.extent1();
  int nrhs          = 1;
  int lda           = m;
  int ldb           = std::max(m, n);
  int lwork         = work.extent0() * work.extent1();
  int info;

  SIERRA_FORTRAN(dgels)(&transC, &m, &n, &nrhs, at.data(), &lda, b.data(), &ldb, work.data(), &lwork, &info);

  // ordinarily A would be passed into lapack, where it would get overwritten.
  // Trash the contents of A so users dont reuse A, in case we change the interface
  // in the future to avoid the transpose
  a.fill(1);

  if (info != 0)
    throw std::runtime_error("dgels returned an error value: " + std::to_string(info));
}

void compute_qr_factorization(Matrix<double>& a, Matrix<double>& work, double* tau)
{
  Matrix<double> at = transpose(a);
  int m             = a.extent0();
  int n             = a.extent1();
  int lda           = m;
  int lwork         = work.extent0() * work.extent1();
  int info;

  SIERRA_FORTRAN(dgeqrf)(&m, &n, at.data(), &lda, tau, work.data(), &lwork, &info);

  std::memcpy(a.data(), at.data(), a.extent0() * a.extent1() * sizeof(double));

  if (info != 0)
    throw std::runtime_error("dgeqrf returned an error value: " + std::to_string(info));
}

void compute_rank_revealing_qr(Matrix<double>& a, Matrix<double>& work, double* tau, int* jpvt)
{
  Matrix<double> at = transpose(a);
  int m             = a.extent0();
  int n             = a.extent1();
  int lda           = m;
  int lwork         = work.extent0() * work.extent1();
  int info;

  SIERRA_FORTRAN(dgeqp3)(&m, &n, at.data(), &lda, jpvt, tau, work.data(), &lwork, &info);

  if (info != 0)
    throw std::runtime_error("dgeqp3 returned an error value: " + std::to_string(info));  

  std::memcpy(a.data(), at.data(), a.extent0() * a.extent1() * sizeof(double));
}

// Given A, the output of dgeqrf (the QR factorization computed by lapack), multiply the
// rhs by Q^T
void multiply_q_transposed(Matrix<double>& a, Matrix<double>& work, double* tau, double* rhs)
{
  char side  = 'L';
  char trans = 'T';
  int m      = a.extent0();
  int n      = 1;
  int k      = m;
  int lda    = m;
  int ldc    = m;
  int lwork  = work.extent0() * work.extent1();
  int info;

  SIERRA_FORTRAN(dormqr)(&side, &trans, &m, &n, &k, a.data(), &lda, tau, rhs, &ldc, work.data(), &lwork, &info);

  if (info != 0)
    throw std::runtime_error("dormqr returned an error value: " + std::to_string(info));
}

// solves an upper triangular system using only the first m rows/columns of the triangular factor
// in A 
void solve_upper_triangular(Matrix<double>& a, double* rhs, int m)
{
  char side    = 'L';
  char uplo    = 'U';
  char transA  = 'N';
  char diag    = 'N';
  //int m        = a.extent0();
  int n        = 1;
  double alpha = 1;
  int lda      = a.extent0();
  int ldb      = a.extent0();

  SIERRA_FORTRAN(dtrsm)(&side, &uplo, &transA, &diag, &m, &n, &alpha, a.data(), &lda, rhs, &ldb);
}

void solve_upper_triangular(Matrix<double>& a, double* rhs)
{
  solve_upper_triangular(a, rhs, a.extent0());
}

void solve_qr_factorization(Matrix<double>& a, Matrix<double>& work, double* tau, double* rhs)
{
  multiply_q_transposed(a, work, tau, rhs);
  solve_upper_triangular(a, rhs);
}

void solve_qrp_factorization(Matrix<double>& a, Matrix<double>& work, double* tau, int* jpvt, double* rhs)
{
  solve_qrp_factorization(a, work, tau, jpvt, rhs, std::min(a.extent0(), a.extent1()));
}

void solve_qrp_factorization(Matrix<double>& a, Matrix<double>& work, double* tau, int* jpvt, double* rhs, int numSingularVals)
{
  multiply_q_transposed(a, work, tau, rhs);
  solve_upper_triangular(a, rhs, numSingularVals);
  for (int i=numSingularVals; i < a.extent1(); ++i)
    rhs[i] = 0;

  assert(work.extent0()*work.extent1() >= a.extent1());
  double* workptr = &(work(0, 0));
  for (int i=0; i < a.extent1(); ++i)
  {
    workptr[jpvt[i]-1] = rhs[i];
  }

  for (int i=0; i < a.extent1(); ++i)
  {
    rhs[i] = workptr[i];
  }
}


void solve_linear_system(Matrix<double>& a, int* ipiv, double* rhs, int nrhs)
{
  int info        = -1;
  int n           = a.extent0();
  const char flag = 'T';

  SIERRA_FORTRAN(dgetrf)(&n, &n, a.data(), &n, ipiv, &info);

  if (info != 0)
    throw std::runtime_error("dgetrf returned an error value: " + std::to_string(info));

  SIERRA_FORTRAN(dgetrs)(&flag, &n, &nrhs, a.data(), &n, ipiv, rhs, &n, &info);

  if (info != 0)
    throw std::runtime_error("dgetrs returned an error value: " + std::to_string(info));
}

void compute_svd_factorization(Matrix<double>& a, double* s, double* u, double* vt, double* work, int lwork)
{
  int m    = a.extent0();
  int n    = a.extent1();
  int lda  = m;
  int ldu  = m;
  int ldvt = n;
  const char jobu = 'A';
  const char jobvt = 'A';
  Matrix<double> at = transpose(a);
  int info;

  SIERRA_FORTRAN(dgesvd)(&jobu, &jobvt, &m, &n, at.data(), &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);

  if (info != 0)
    throw std::runtime_error("dgesvd returned an error value: " + std::to_string(info));

  a.fill(-1);  // if a had been passed into lapack, it would have been overwritten with garbage.
               // Do that here so users dont start expecting a to be unchanged
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
