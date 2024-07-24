// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_IMPL_SYSTEMBLAS_HPP
#define TSQR_IMPL_SYSTEMBLAS_HPP

#include "Tsqr_ConfigDefs.hpp"
#include "Tsqr_Impl_RawBlas.hpp"
#include "Teuchos_BLAS_types.hpp"
#include <complex>

namespace TSQR {
namespace Impl {

template<class Scalar>
class SystemBlas {};

#define TSQR_IMPL_SYSTEMBLAS_DECL( Scalar ) \
template<> \
class SystemBlas<Scalar> : public RawBlas<Scalar> { \
public: \
  using value_type = Scalar; \
  \
  ~SystemBlas() = default; \
  \
  void \
  matrix_matrix_product(const char transa, const char transb, \
                        const int m, const int n, const int k, \
                        const value_type& alpha, \
                        const value_type A[], const int lda, \
                        const value_type B[], const int ldb, \
                        const value_type& beta, \
                        value_type C[], const int ldc) const override; \
  \
  void \
  GEMM(const Teuchos::ETransp transa, const Teuchos::ETransp transb, \
       const int m, const int n, const int k,                        \
       const value_type& alpha,                                      \
       const value_type A[], const int lda,                          \
       const value_type B[], const int ldb,                          \
       const value_type& beta,                                         \
       value_type C[], const int ldc) const;                  \
  \
  void \
  triangular_matrix_matrix_solve(const char side, const char uplo, \
                                 const char transa, const char diag, \
                                 const int m, const int n, \
                                 const value_type& alpha, \
                                 const value_type A[], const int lda, \
                                 value_type B[], const int ldb) const override; \
  \
  void \
  TRSM(const Teuchos::ESide side, const Teuchos::EUplo uplo, \
       const Teuchos::ETransp transa, const Teuchos::EDiag diag, \
       const int m, const int n, \
       const value_type& alpha, \
       const value_type A[], const int lda, \
       value_type B[], const int ldb) const; \
};

TSQR_IMPL_SYSTEMBLAS_DECL( float )
TSQR_IMPL_SYSTEMBLAS_DECL( double )

#ifdef HAVE_TPETRATSQR_COMPLEX
TSQR_IMPL_SYSTEMBLAS_DECL( std::complex<float> )
TSQR_IMPL_SYSTEMBLAS_DECL( std::complex<double> )
#endif // HAVE_TPETRATSQR_COMPLEX

} // namespace Impl
} // namespace TSQR

#endif // TSQR_IMPL_SYSTEMBLAS_HPP
