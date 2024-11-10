// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_IMPL_RAWBLAS_HPP
#define TSQR_IMPL_RAWBLAS_HPP

namespace TSQR {
namespace Impl {

/// \brief "Raw" local BLAS interface.
///
/// Subclass and specialize this interface as needed.
///
/// The methods are instance methods so that subclass instances may
/// have state.  For example, a cuBLAS implementation would have a
/// CUDA stream instance (cudaStream_t) and some kind of handle.
///
/// Methods are virtual because they are meant to be called from host,
/// even if they run on device with pointers to device data.
template<class Scalar>
class RawBlas {
public:
  using value_type = Scalar;

  virtual ~RawBlas() = default;

  //! Corresponds to BLAS _GEMM.
  virtual void
  matrix_matrix_product(const char transa, const char transb,
                        const int m, const int n, const int k,
                        const value_type& alpha,
                        const value_type A[], const int lda,
                        const value_type B[], const int ldb,
                        const value_type& beta,
                        value_type C[], const int ldc) const = 0;

  //! Corresponds to BLAS _TRSM.
  virtual void
  triangular_matrix_matrix_solve(const char side, const char uplo,
                                 const char transa, const char diag,
                                 const int m, const int n,
                                 const value_type& alpha,
                                 const value_type A[], const int lda,
                                 value_type B[], const int ldb) const = 0;
};

} // namespace Impl
} // namespace TSQR

#endif // TSQR_IMPL_RAWBLAS_HPP
