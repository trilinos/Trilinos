// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_IMPL_CUSOLVER_HPP
#define TSQR_IMPL_CUSOLVER_HPP

#include "TpetraTSQR_config.h"
#if defined(HAVE_TPETRATSQR_CUBLAS) && defined(HAVE_TPETRATSQR_CUSOLVER)
#include "Tsqr_Impl_CuSolverHandle_fwd.hpp"
#if defined(HAVE_TPETRATSQR_COMPLEX)
#  include <complex>
#endif // HAVE_TPETRATSQR_COMPLEX
#include "Tsqr_Impl_RawQR.hpp"

namespace TSQR {
namespace Impl {

template<class Scalar>
class CuSolver : public RawQR<Scalar> {
public:
  CuSolver(int* const info); // use default cuSOLVER handle

  CuSolver(const std::shared_ptr<CuSolverHandle>& handle,
           int* const info);

  virtual bool wants_device_memory () const { return true; }

  int
  compute_QR_lwork(const int nrows,
                   const int ncols,
                   Scalar A_raw[],
                   const int lda) const override;

  void
  compute_QR(const int nrows,
             const int ncols,
             Scalar A[],
             const int lda,
             Scalar tau[],
             Scalar work[],
             const int lwork) const override;

  int
  apply_Q_factor_lwork(const char side,
                       const char trans,
                       const int nrows,
                       const int ncols_C,
                       const int ncols_Q,
                       const Scalar Q[],
                       const int ldq,
                       const Scalar tau[],
                       Scalar C[],
                       const int ldc) const override;

  void
  apply_Q_factor(const char side,
                 const char trans,
                 const int nrows,
                 const int ncols_C,
                 const int ncols_Q,
                 const Scalar Q[],
                 const int ldq,
                 const Scalar tau[],
                 Scalar C[],
                 const int ldc,
                 Scalar work[],
                 const int lwork) const override;

  int
  compute_explicit_Q_lwork(const int m, const int n, const int k,
                           Scalar A[], const int lda,
                           const Scalar tau[]) const override;

  void
  compute_explicit_Q(const int m, const int n, const int k,
                     Scalar A[], const int lda,
                     const Scalar tau[],
                     Scalar work[], const int lwork) const override;

private:
  std::shared_ptr<CuSolverHandle> handle_;
  int* info_; // DEVICE MEMORY
};

extern template class CuSolver<double>;
extern template class CuSolver<float>;
#if defined(HAVE_TPETRATSQR_COMPLEX)
extern template class CuSolver<std::complex<double>>;
extern template class CuSolver<std::complex<float>>;
#endif // defined(HAVE_TPETRATSQR_COMPLEX)

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS && HAVE_TPETRATSQR_CUSOLVER
#endif // TSQR_IMPL_CUSOLVER_HPP
