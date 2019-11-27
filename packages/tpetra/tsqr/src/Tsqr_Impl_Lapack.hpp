#ifndef TSQR_IMPL_LAPACK_HPP
#define TSQR_IMPL_LAPACK_HPP

#include "Tsqr_ConfigDefs.hpp"
#include "Tsqr_Impl_RawQR.hpp"
#include <complex>

namespace TSQR {
namespace Impl {

template<class Scalar>
class Lapack {};

// CombineNative needs LARFG, but it's not properly part of RawQR.
// RawQR needs to be able to wrap lots of different functions,
// including whatever cuSOLVER provides.  It doesn't make sense to
// launch a device kernel from host for ever column of the matrix,
// especially not when cuSOLVER already has all the needed QR
// factorization and apply Q factor functions.

#define TSQR_IMPL_LAPACK_DECL( Scalar ) \
template<> \
class Lapack<Scalar> : public RawQR<Scalar> { \
public: \
  using value_type = Scalar; \
  using magnitude_type = decltype(std::abs(Scalar{})); \
  \
  ~Lapack() = default; \
  \
  void \
  compute_QR(const int m, const int n, value_type A[], \
             const int lda, value_type TAU[], value_type WORK[], \
             const int lwork) const override; \
  \
  void \
  apply_Q_factor(const char SIDE, const char TRANS, \
                 const int m, const int n, const int k, \
                 const value_type A[], const int lda, \
                 const value_type TAU[], \
                 value_type C[], const int ldc, \
                 value_type WORK[], const int lwork) const override; \
  \
  void \
  compute_explicit_Q(const int m, const int n, const int k, \
                     value_type A[], const int lda, \
                     const value_type TAU[], value_type WORK[], \
                     const int lwork) const override; \
  \
  void \
  GESVD(const char JOBU, const char JOBVT, \
        const int m, const int n, \
        value_type A[], const int lda, \
        magnitude_type S[], value_type U[], const int ldu, \
        value_type V[], const int ldv, \
        value_type WORK[], const int lwork, \
        magnitude_type RWORK[]) const; \
  \
  void \
  LARFG(const int n, value_type& alpha, value_type x[], \
        const int incx, value_type& tau) const; \
  \
  void \
  POTRF(const char UPLO, const int n, \
        value_type A[], const int lda) const; \
  \
  void \
  LARNV(const int idist, int seed[], const int n, \
        value_type v[]) const; \
};

TSQR_IMPL_LAPACK_DECL( float )
TSQR_IMPL_LAPACK_DECL( double )

#ifdef HAVE_KOKKOSTSQR_COMPLEX
TSQR_IMPL_LAPACK_DECL( std::complex<float> )
TSQR_IMPL_LAPACK_DECL( std::complex<double> )
#endif // HAVE_KOKKOSTSQR_COMPLEX

} // namespace Impl
} // namespace TSQR

#endif // TSQR_IMPL_LAPACK_HPP
