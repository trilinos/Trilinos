#ifndef __TACHO_BLAS_EXTERNAL_HPP__
#define __TACHO_BLAS_EXTERNAL_HPP__

/// \file  Tacho_Blas_External.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Core.hpp" // CUDA specialization

#if defined(KOKKOS_ENABLE_CUDA) 
#include "cublas_v2.h"
#endif

namespace Tacho {

  template<typename T>
  struct Blas {
    static 
    int gemv(const char trans, 
             int m, int n, 
             const T alpha, 
             const T *a, int lda,
             const T *b, int ldb,
             const T beta,
             /* */ T *c, int ldc);
#if defined (KOKKOS_ENABLE_CUDA)
    static 
    int gemv(cublasHandle_t handle,
             const cublasOperation_t trans,
             int m, int n, 
             const T alpha, 
             const T *a, int lda,
             const T *b, int ldb,
             const T beta,
             /* */ T *c, int ldc);
#endif

    static 
    int trsv(const char uplo, const char transa, const char diag, 
             int m, 
             const T *a, int lda,
             /* */ T *b, int ldb);
#if defined (KOKKOS_ENABLE_CUDA)
    static 
    int trsv(cublasHandle_t handle,
             const cublasFillMode_t uplo, const cublasOperation_t transa, const cublasDiagType_t diag,
             int m, 
             const T *a, int lda,
             /* */ T *b, int ldb);
#endif

    static 
    int gemm(const char transa, const char transb, 
             int m, int n, int k,
             const T alpha, 
             const T *a, int lda,
             const T *b, int ldb,
             const T beta,
             /* */ T *c, int ldc);

    static 
    int herk(const char uplo, const char trans, 
             int n, int k,
             const T alpha, 
             const T *a, int lda,
             const T beta,
             /* */ T *c, int ldc);
#if defined (KOKKOS_ENABLE_CUDA)
    static 
    int herk(cublasHandle_t handle,
             const cublasFillMode_t uplo, const cublasOperation_t trans,
             int n, int k,
             const T alpha, 
             const T *a, int lda,
             const T beta,
             /* */ T *c, int ldc);
#endif

    static 
    int trsm(const char side, const char uplo, const char transa, const char diag,
             int m, int n, 
             const T alpha, 
             const T *a, int lda,
             /* */ T *b, int ldb);
#if defined (KOKKOS_ENABLE_CUDA)
    static 
    int trsm(cublasHandle_t handle,
             const cublasSideMode_t side, const cublasFillMode_t uplo,
             const cublasOperation_t trans, const cublasDiagType_t diag,
             int m, int n, 
             const T alpha, 
             const T *a, int lda,
             /* */ T *b, int ldb);
#endif
  };

}

#endif
