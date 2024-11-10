// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_Impl_SystemBlas.hpp"
#include "Teuchos_BLAS.hpp"

namespace TSQR {
namespace Impl {

#define TSQR_IMPL_SYSTEMBLAS_IMPL( Scalar ) \
void SystemBlas<Scalar>:: \
matrix_matrix_product(const char transa, const char transb, \
                      const int m, const int n, const int k, \
                      const value_type& alpha, \
                      const value_type A[], const int lda, \
                      const value_type B[], const int ldb, \
                      const value_type& beta, \
                      value_type C[], const int ldc) const \
{ \
  const Teuchos::ETransp transa_enum = \
    (transa == 'C' || transa == 'c') ? \
    Teuchos::CONJ_TRANS : \
    ((transa == 'T' || transa == 't') ? \
     Teuchos::TRANS : \
     Teuchos::NO_TRANS); \
  const Teuchos::ETransp transb_enum = \
    (transb == 'C' || transb == 'c') ? \
    Teuchos::CONJ_TRANS : \
    ((transb == 'T' || transb == 't') ? \
     Teuchos::TRANS : \
     Teuchos::NO_TRANS); \
  GEMM(transa_enum, transb_enum, m, n, k, \
       alpha, A, lda, B, ldb, beta, C, ldc); \
} \
  \
void SystemBlas<Scalar>:: \
GEMM(const Teuchos::ETransp transa, const Teuchos::ETransp transb, \
     const int m, const int n, const int k, \
     const value_type& alpha, const value_type A[], const int lda, \
     const value_type B[], const int ldb, \
     const value_type& beta, value_type C[], const int ldc) const \
{ \
  Teuchos::BLAS<int, value_type> blas; \
  blas.GEMM(transa, transb, m, n, k, \
            alpha, A, lda, B, ldb, beta, C, ldc); \
} \
  \
void SystemBlas<Scalar>:: \
triangular_matrix_matrix_solve(const char side, const char uplo, \
                               const char transa, const char diag, \
                               const int m, const int n, \
                               const value_type& alpha, \
                               const value_type A[], const int lda, \
                               value_type B[], const int ldb) const \
{ \
  const Teuchos::ESide side_enum = \
    (side == 'L' || side == 'l') ? \
    Teuchos::LEFT_SIDE : \
    Teuchos::RIGHT_SIDE; \
  const Teuchos::EUplo uplo_enum = \
    (uplo == 'U' || uplo == 'u') ? \
    Teuchos::UPPER_TRI : \
    ((uplo == 'L' || uplo == 'l') ? \
     Teuchos::LOWER_TRI : \
     Teuchos::UNDEF_TRI); \
  const Teuchos::ETransp transa_enum = \
    (transa == 'C' || transa == 'c') ? \
    Teuchos::CONJ_TRANS : \
    ((transa == 'T' || transa == 't') ? \
     Teuchos::TRANS : \
     Teuchos::NO_TRANS); \
  const Teuchos::EDiag diag_enum = \
    (diag == 'U' || diag == 'u') ? \
    Teuchos::UNIT_DIAG : \
    Teuchos::NON_UNIT_DIAG; \
  TRSM(side_enum, uplo_enum, transa_enum, diag_enum, \
       m, n, alpha, A, lda, B, ldb); \
} \
  \
void SystemBlas<Scalar>:: \
TRSM(const Teuchos::ESide side, const Teuchos::EUplo uplo, \
     const Teuchos::ETransp transa, const Teuchos::EDiag diag, \
     const int m, const int n, \
     const value_type& alpha, \
     const value_type A[], const int lda, \
     value_type B[], const int ldb) const \
{ \
  Teuchos::BLAS<int, value_type> blas; \
  blas.TRSM(side, uplo, transa, diag, \
            m, n, alpha, A, lda, B, ldb); \
}

TSQR_IMPL_SYSTEMBLAS_IMPL( float )
TSQR_IMPL_SYSTEMBLAS_IMPL( double )

#ifdef HAVE_TPETRATSQR_COMPLEX
TSQR_IMPL_SYSTEMBLAS_IMPL( std::complex<float> )
TSQR_IMPL_SYSTEMBLAS_IMPL( std::complex<double> )
#endif // HAVE_TPETRATSQR_COMPLEX

} // namespace Impl
} // namespace TSQR
