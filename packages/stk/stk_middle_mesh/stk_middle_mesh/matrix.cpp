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

// TODO: LAPACKE_ROW_MAJOR typically causes Lapack to malloc temporary arrays and
//       transpose things
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

void solve_upper_triangular(Matrix<double>& a, double* rhs)
{
  char side    = 'L';
  char uplo    = 'U';
  char transA  = 'N';
  char diag    = 'N';
  int m        = a.extent0();
  int n        = 1;
  double alpha = 1;
  int lda      = m;
  int ldb      = m;

  SIERRA_FORTRAN(dtrsm)(&side, &uplo, &transA, &diag, &m, &n, &alpha, a.data(), &lda, rhs, &ldb);
}

void solve_qr_factorization(Matrix<double>& a, Matrix<double>& work, double* tau, double* rhs)
{
  multiply_q_transposed(a, work, tau, rhs);
  solve_upper_triangular(a, rhs);
}

void solve_linear_system(Matrix<double>& a, int* ipiv, double* rhs)
{
  int info        = -1;
  int val         = 1;
  int n           = a.extent0();
  const char flag = 'T';

  SIERRA_FORTRAN(dgetrf)(&n, &n, a.data(), &n, ipiv, &info);

  if (info != 0)
    throw std::runtime_error("dgetrf returned an error value: " + std::to_string(info));

  SIERRA_FORTRAN(dgetrs)(&flag, &n, &val, a.data(), &n, ipiv, rhs, &n, &info);

  if (info != 0)
    throw std::runtime_error("dgetrs returned an error value: " + std::to_string(info));
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
