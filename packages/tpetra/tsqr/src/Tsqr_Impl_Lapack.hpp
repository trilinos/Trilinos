// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_IMPL_LAPACK_HPP
#define TSQR_IMPL_LAPACK_HPP

#include "Tsqr_ConfigDefs.hpp"
#include "Tsqr_Impl_RawQR.hpp"
#include <complex>

namespace TSQR {
namespace Impl {

/// \brief Implementation of RawQR that uses the system's LAPACK
///   library via Teuchos::LAPACK.
///
/// This class provides functions not in RawQR for the sake of
/// CombineNative.  CombineNative needs LARFG, but it's not properly
/// part of RawQR.  It doesn't make sense to launch a device kernel
/// from host for every column of the matrix, especially not when
/// cuSOLVER already has all the needed QR factorization and apply Q
/// factor functions.
template<class Scalar>
class Lapack : public RawQR<Scalar> {
public:
  using value_type = Scalar;
  using magnitude_type = decltype(std::abs(Scalar{}));

  // NOTE (mfh 22 Dec 2019) I would normally write "= default;" here,
  // but Intel 17 appears to have a bug that requires an explicit
  // nondefault definition.  See discussion here:
  //
  // https://github.com/trilinos/Trilinos/pull/6488#issuecomment-568351758
  ~Lapack() override {}

  int
  compute_QR_lwork(const int m, const int n,
                   value_type A[], const int lda) const override;

  void
  compute_QR(const int m, const int n, value_type A[],
             const int lda, value_type TAU[], value_type WORK[],
             const int lwork) const override;

  int
  apply_Q_factor_lwork(const char SIDE, const char TRANS,
                       const int m, const int n, const int k,
                       const value_type A[], const int lda,
                       const value_type TAU[],
                       value_type C[], const int ldc) const override;

  void
  apply_Q_factor(const char SIDE, const char TRANS,
                 const int m, const int n, const int k,
                 const value_type A[], const int lda,
                 const value_type TAU[],
                 value_type C[], const int ldc,
                 value_type WORK[], const int lwork) const override;

  int
  compute_explicit_Q_lwork(const int m, const int n, const int k,
                           value_type A[], const int lda,
                           const value_type TAU[]) const override;

  void
  compute_explicit_Q(const int m, const int n, const int k,
                     value_type A[], const int lda,
                     const value_type TAU[], value_type WORK[],
                     const int lwork) const override;

  void
  GESVD(const char JOBU, const char JOBVT,
        const int m, const int n,
        value_type A[], const int lda,
        magnitude_type S[], value_type U[], const int ldu,
        value_type V[], const int ldv,
        value_type WORK[], const int lwork,
        magnitude_type RWORK[]) const;

  void
  LARFG(const int n, value_type& alpha, value_type x[],
        const int incx, value_type& tau) const;

  void
  POTRF(const char UPLO, const int n,
        value_type A[], const int lda) const;

  void
  LARNV(const int idist, int seed[], const int n,
        value_type v[]) const;
};

extern template class Lapack<float>;
extern template class Lapack<double>;

#ifdef HAVE_TPETRATSQR_COMPLEX
extern template class Lapack<std::complex<float>>;
extern template class Lapack<std::complex<double>>;
#endif // HAVE_TPETRATSQR_COMPLEX

} // namespace Impl
} // namespace TSQR

#endif // TSQR_IMPL_LAPACK_HPP
