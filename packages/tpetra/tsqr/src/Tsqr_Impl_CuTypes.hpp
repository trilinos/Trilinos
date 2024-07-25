// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_IMPL_CUTYPES_HPP
#define TSQR_IMPL_CUTYPES_HPP

#include "TpetraTSQR_config.h"
#if defined(HAVE_TPETRATSQR_CUBLAS) && defined(HAVE_TPETRATSQR_CUSOLVER)
#include <cublas_v2.h> // for cublasSideMode_t etc.
#include <cusolverDn.h>
#if defined(HAVE_TPETRATSQR_COMPLEX)
#  include <complex>
#endif // HAVE_TPETRATSQR_COMPLEX

namespace TSQR {
namespace Impl {

template<class Scalar>
struct CudaValue {};

template<>
struct CudaValue<double> {
  using type = double;

  static type makeValue (const double x) {
    return x;
  }

  static bool arrayCorrectlyAligned (const double* const /* x */) {
    return true;
  }
};

template<>
struct CudaValue<float> {
  using type = float;

  static type makeValue (const float x) {
    return x;
  }

  static bool arrayCorrectlyAligned (const double* const /* x */) {
    return true;
  }
};

#if defined(HAVE_TPETRATSQR_COMPLEX)
// FIXME (mfh 10 Dec 2019) CUDA's built-in complex types must be
// aligned to the whole type, not just to double or float (as with
// std::complex or (currently) Kokkos::complex).
template<>
struct CudaValue<std::complex<double>> {
  using type = cuDoubleComplex;

  static type makeValue (const std::complex<double> x) {
    return make_cuDoubleComplex (std::real (x), std::imag (x));
  }

  static bool
  arrayCorrectlyAligned (const std::complex<double>* const x)
  {
    // CUDA requires arrays of complex to be aligned to the full type,
    // not just to one of the two numbers (as with std::complex).
    constexpr size_t requiredAlignment =
      sizeof (std::complex<double>);
    return x == nullptr ||
      reinterpret_cast<size_t> (x) % requiredAlignment == 0;
  }
};

template<>
struct CudaValue<std::complex<float>> {
  using type = cuFloatComplex;

  static type makeValue (const std::complex<float> x) {
    return make_cuFloatComplex (std::real (x), std::imag (x));
  }

  static bool
  arrayCorrectlyAligned (const std::complex<float>* const x)
  {
    // CUDA requires arrays of complex to be aligned to the full type,
    // not just to one of the two numbers (as with std::complex).
    constexpr size_t requiredAlignment =
      sizeof (std::complex<float>);
    return x == nullptr ||
      reinterpret_cast<size_t> (x) % requiredAlignment == 0;
  }
};
#endif // defined(HAVE_TPETRATSQR_COMPLEX)

cublasSideMode_t cuBlasSide (const char side);
cublasOperation_t cuBlasTrans (const char trans);

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS && HAVE_TPETRATSQR_CUSOLVER
#endif // TSQR_IMPL_CUTYPES_HPP
