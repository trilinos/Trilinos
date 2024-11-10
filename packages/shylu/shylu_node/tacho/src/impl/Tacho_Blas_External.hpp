// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_BLAS_EXTERNAL_HPP__
#define __TACHO_BLAS_EXTERNAL_HPP__

/// \file  Tacho_Blas_External.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Core.hpp"

#if defined(KOKKOS_ENABLE_CUDA)
#define TACHO_ENABLE_CUBLAS
#endif

#if defined(KOKKOS_ENABLE_HIP)
#define TACHO_ENABLE_ROCBLAS
#endif

#if defined(TACHO_ENABLE_CUBLAS)
#include "cublas_v2.h"
#endif

#if defined(TACHO_ENABLE_ROCBLAS)
#include "rocblas/rocblas.h"
#endif

namespace Tacho {

template <typename T> struct Blas {
  static int gemv(const char trans, int m, int n, const T alpha, const T *a, int lda, const T *b, int ldb, const T beta,
                  /* */ T *c, int ldc);
#if defined(TACHO_ENABLE_CUBLAS)
  static int gemv(cublasHandle_t handle, const cublasOperation_t trans, int m, int n, const T alpha, const T *a,
                  int lda, const T *b, int ldb, const T beta,
                  /* */ T *c, int ldc);
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
  static int gemv(rocblas_handle handle, const rocblas_operation trans, int m, int n, const T alpha, const T *a,
                  int lda, const T *b, int ldb, const T beta,
                  /* */ T *c, int ldc);
#endif

  static int trsv(const char uplo, const char transa, const char diag, int m, const T *a, int lda,
                  /* */ T *b, int ldb);
#if defined(TACHO_ENABLE_CUBLAS)
  static int trsv(cublasHandle_t handle, const cublasFillMode_t uplo, const cublasOperation_t transa,
                  const cublasDiagType_t diag, int m, const T *a, int lda,
                  /* */ T *b, int ldb);
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
  static int trsv(rocblas_handle handle, const rocblas_fill uplo, const rocblas_operation transa,
                  const rocblas_diagonal diag, int m, const T *a, int lda,
                  /* */ T *b, int ldb);
#endif

  static int gemm(const char transa, const char transb, int m, int n, int k, const T alpha, const T *a, int lda,
                  const T *b, int ldb, const T beta,
                  /* */ T *c, int ldc);
#if defined(TACHO_ENABLE_CUBLAS)
  static int gemm(cublasHandle_t handle, const cublasOperation_t transa, const cublasOperation_t transb, int m, int n,
                  int k, const T alpha, const T *a, int lda, const T *b, int ldb, const T beta,
                  /* */ T *c, int ldc);
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
  static int gemm(rocblas_handle handle, const rocblas_operation transa, const rocblas_operation transb, int m, int n,
                  int k, const T alpha, const T *a, int lda, const T *b, int ldb, const T beta,
                  /* */ T *c, int ldc);
#endif

  static int herk(const char uplo, const char trans, int n, int k, const T alpha, const T *a, int lda, const T beta,
                  /* */ T *c, int ldc);
#if defined(TACHO_ENABLE_CUBLAS)
  static int herk(cublasHandle_t handle, const cublasFillMode_t uplo, const cublasOperation_t trans, int n, int k,
                  const T alpha, const T *a, int lda, const T beta,
                  /* */ T *c, int ldc);
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
  static int herk(rocblas_handle handle, const rocblas_fill uplo, const rocblas_operation trans, int n, int k,
                  const T alpha, const T *a, int lda, const T beta,
                  /* */ T *c, int ldc);
#endif

  static int trsm(const char side, const char uplo, const char transa, const char diag, int m, int n, const T alpha,
                  const T *a, int lda,
                  /* */ T *b, int ldb);
#if defined(TACHO_ENABLE_CUBLAS)
  static int trsm(cublasHandle_t handle, const cublasSideMode_t side, const cublasFillMode_t uplo,
                  const cublasOperation_t trans, const cublasDiagType_t diag, int m, int n, const T alpha, const T *a,
                  int lda,
                  /* */ T *b, int ldb);
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
  static int trsm(rocblas_handle handle, const rocblas_side side, const rocblas_fill uplo,
                  const rocblas_operation trans, const rocblas_diagonal diag, int m, int n, const T alpha, const T *a,
                  int lda,
                  /* */ T *b, int ldb);
#endif
};

} // namespace Tacho

#endif
