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
struct CuSolverValue {};

template<>
struct CuSolverValue<double> {
  using type = double;
};

template<>
struct CuSolverValue<float> {
  using type = float;
};

#if defined(HAVE_TPETRATSQR_COMPLEX)
// FIXME (mfh 10 Dec 2019) CUDA's built-in complex types must be
// aligned to the whole type, not just to double or float (as with
// std::complex or (currently) Kokkos::complex).
template<>
struct CuSolverValue<std::complex<double>> {
  using type = cuDoubleComplex;
};

template<>
struct CuSolverValue<std::complex<float>> {
  using type = cuComplex;
};
#endif // defined(HAVE_TPETRATSQR_COMPLEX)

cublasSideMode_t cuBlasSide (const char side);
cublasOperation_t cuBlasTrans (const char trans);

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS && HAVE_TPETRATSQR_CUSOLVER
#endif // TSQR_IMPL_CUTYPES_HPP
