// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_IMPL_CUBLAS_HPP
#define TSQR_IMPL_CUBLAS_HPP

#include "TpetraTSQR_config.h"
#if defined(HAVE_TPETRATSQR_CUBLAS)
#  include "Tsqr_Impl_CuBlasHandle_fwd.hpp"
#  if defined(HAVE_TPETRATSQR_COMPLEX)
#    include <complex>
#  endif // HAVE_TPETRATSQR_COMPLEX

namespace TSQR {
namespace Impl {

template<class Scalar>
class CuBlas {
public:
  // Use the default handle.
  CuBlas ();
  CuBlas (const std::shared_ptr<CuBlasHandle>& handle);

  void
  gemm (const char transa,
        const char transb,
        const int m, const int n, const int k,
        const Scalar alpha,
        const Scalar* A, const int lda,
        const Scalar* B, const int ldb,
        const Scalar beta,
        Scalar* C, const int ldc);

private:
  std::shared_ptr<CuBlasHandle> handle_;
};

extern template class CuBlas<double>;
extern template class CuBlas<float>;
#if defined(HAVE_TPETRATSQR_COMPLEX)
extern template class CuBlas<std::complex<double>>;
extern template class CuBlas<std::complex<float>>;
#endif // defined(HAVE_TPETRATSQR_COMPLEX)

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS
#endif // TSQR_IMPL_CUBLAS_HPP
