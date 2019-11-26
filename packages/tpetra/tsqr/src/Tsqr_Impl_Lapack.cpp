#include "Tsqr_Impl_Lapack.hpp"
#include "Teuchos_LAPACK.hpp"
#include <sstream>
#include <stdexcept>

namespace TSQR {
namespace Impl {

#define TSQR_IMPL_LAPACK_IMPL( Scalar ) \
void Lapack<Scalar>:: \
LARNV(const int idist, int seed[], const int n, \
      value_type v[]) const \
{ \
  Teuchos::LAPACK<int, value_type> lapack; \
  lapack.LARNV(idist, seed, n, v); \
} \
  \
void Lapack<Scalar>:: \
POTRF(const char UPLO, const int n, \
      value_type A[], const int lda) const \
{ \
  Teuchos::LAPACK<int, value_type> lapack; \
  int info = 0; \
  lapack.POTRF(UPLO, n, A, lda, &info); \
  if (info != 0) { \
    std::ostringstream os; \
    os << "LAPACK POTRF (Cholesky factorization) " \
       << "failed with INFO = " << info << "."; \
    throw std::logic_error (os.str ()); \
  } \
} \
  \
void Lapack<Scalar>:: \
GESVD(const char JOBU, const char JOBVT, \
      const int m, const int n, \
      value_type A[], const int lda, \
      magnitude_type S[], value_type U[], const int ldu, \
      value_type V[], const int ldv, \
      value_type WORK[], const int lwork, \
      magnitude_type RWORK[]) const \
{ \
  Teuchos::LAPACK<int, value_type> lapack; \
  int info = 0; \
  lapack.GESVD(JOBU, JOBVT, m, n, A, lda, S, \
               U, ldu, V, ldv, WORK, lwork, RWORK, &info); \
  if (info != 0) { \
    std::ostringstream os; \
    os << "LAPACK GESVD (singular value decomposition) " \
       << "failed with INFO = " << info << "."; \
    throw std::logic_error (os.str ()); \
  } \
} \
  \
void Lapack<Scalar>:: \
LARFG(const int n, value_type& alpha, value_type x[], \
      const int incx, value_type& tau) const \
{ \
  Teuchos::LAPACK<int, value_type> lapack; \
  lapack.LARFG(n, &alpha, x, incx, &tau); \
} \
  \
void Lapack<Scalar>:: \
compute_QR(const int m, const int n, value_type A[], const int lda, \
           value_type TAU[], value_type WORK[], const int lwork) const \
{ \
  Teuchos::LAPACK<int, value_type> lapack; \
  int info = 0; \
  lapack.GEQRF(m, n, A, lda, TAU, WORK, lwork, &info); \
  if (info != 0) { \
    std::ostringstream os; \
    os << "LAPACK GEQRF (QR factorization) failed with INFO = " \
       << info << "."; \
    throw std::logic_error (os.str()); \
  } \
} \
  \
void Lapack<Scalar>:: \
apply_Q_factor(const char SIDE, const char TRANS, \
               const int m, const int n, const int k, \
               const value_type A[], const int lda, \
               const value_type TAU[], \
               value_type C[], const int ldc, \
               value_type WORK[], const int lwork) const \
{ \
  Teuchos::LAPACK<int, value_type> lapack; \
  int info = 0; \
  value_type* A_nc = const_cast<value_type*>(A); \
  lapack.UNMQR(SIDE, TRANS, m, n, k, A_nc, lda, TAU, C, ldc, WORK, \
               lwork, &info); \
  if (info != 0) { \
    std::ostringstream os; \
    os << "LAPACK UNMQR (apply Q factor from GEQRF) failed with " \
       "INFO = " << info << "."; \
    throw std::logic_error (os.str()); \
  } \
} \
  \
void Lapack<Scalar>:: \
compute_explicit_Q(const int m, const int n, const int k, \
                   value_type A[], const int lda, \
                   const value_type TAU[], value_type WORK[], \
                   const int lwork) const \
{ \
  Teuchos::LAPACK<int, value_type> lapack; \
  int info = 0; \
  lapack.UNGQR(m, n, k, A, lda, TAU, WORK, lwork, &info); \
  if (info != 0) { \
    std::ostringstream os; \
    os << "LAPACK UNGQR (compute explicit Q factor from GEQRF) " \
      "failed with INFO = " << info << "."; \
    throw std::logic_error (os.str()); \
  } \
}

TSQR_IMPL_LAPACK_IMPL( float )
TSQR_IMPL_LAPACK_IMPL( double )

#ifdef HAVE_KOKKOSTSQR_COMPLEX
TSQR_IMPL_LAPACK_IMPL( std::complex<float> )
TSQR_IMPL_LAPACK_IMPL( std::complex<double> )
#endif // HAVE_KOKKOSTSQR_COMPLEX

} // namespace Impl
} // namespace TSQR
