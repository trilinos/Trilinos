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
/// \file  Tacho_Blas_External.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_External.hpp"
#include "Tacho_config.h"

#include "Kokkos_Core.hpp"

extern "C" {
///
/// Gemv
///

void F77_BLAS_MANGLE(sgemv, SGEMV)(const char *, int *, int *, const float *, const float *, int *, const float *,
                                   int *, const float *,
                                   /* */ float *, int *);
void F77_BLAS_MANGLE(dgemv, DGEMV)(const char *, int *, int *, const double *, const double *, int *, const double *,
                                   int *, const double *,
                                   /* */ double *, int *);
void F77_BLAS_MANGLE(cgemv, CGEMV)(const char *, int *, int *, const Kokkos::complex<float> *,
                                   const Kokkos::complex<float> *, int *, const Kokkos::complex<float> *, int *,
                                   const Kokkos::complex<float> *,
                                   /* */ Kokkos::complex<float> *, int *);
void F77_BLAS_MANGLE(zgemv, ZGEMV)(const char *, int *, int *, const Kokkos::complex<double> *,
                                   const Kokkos::complex<double> *, int *, const Kokkos::complex<double> *, int *,
                                   const Kokkos::complex<double> *,
                                   /* */ Kokkos::complex<double> *, int *);

///
/// Trsv
///

void F77_BLAS_MANGLE(strsv, STRSV)(const char *, const char *, const char *, int *, const float *, int *,
                                   /* */ float *, int *);
void F77_BLAS_MANGLE(dtrsv, DTRSV)(const char *, const char *, const char *, int *, const double *, int *,
                                   /* */ double *, int *);
void F77_BLAS_MANGLE(ctrsv, CTRSV)(const char *, const char *, const char *, int *, const Kokkos::complex<float> *,
                                   int *,
                                   /* */ Kokkos::complex<float> *, int *);
void F77_BLAS_MANGLE(ztrsv, ZTRSV)(const char *, const char *, const char *, int *, const Kokkos::complex<double> *,
                                   int *,
                                   /* */ Kokkos::complex<double> *, int *);

///
/// Gemm
///

void F77_BLAS_MANGLE(sgemm, SGEMM)(const char *, const char *, int *, int *, int *, const float *, const float *, int *,
                                   const float *, int *, const float *,
                                   /* */ float *, int *);
void F77_BLAS_MANGLE(dgemm, DGEMM)(const char *, const char *, int *, int *, int *, const double *, const double *,
                                   int *, const double *, int *, const double *,
                                   /* */ double *, int *);
void F77_BLAS_MANGLE(cgemm, CGEMM)(const char *, const char *, int *, int *, int *, const Kokkos::complex<float> *,
                                   const Kokkos::complex<float> *, int *, const Kokkos::complex<float> *, int *,
                                   const Kokkos::complex<float> *,
                                   /* */ Kokkos::complex<float> *, int *);
void F77_BLAS_MANGLE(zgemm, ZGEMM)(const char *, const char *, int *, int *, int *, const Kokkos::complex<double> *,
                                   const Kokkos::complex<double> *, int *, const Kokkos::complex<double> *, int *,
                                   const Kokkos::complex<double> *,
                                   /* */ Kokkos::complex<double> *, int *);

///
/// Herk
///

void F77_BLAS_MANGLE(ssyrk, SSYRK)(const char *, const char *, int *, int *, const float *, const float *, int *,
                                   const float *,
                                   /* */ float *, int *);
void F77_BLAS_MANGLE(dsyrk, DSYRK)(const char *, const char *, int *, int *, const double *, const double *, int *,
                                   const double *,
                                   /* */ double *, int *);
void F77_BLAS_MANGLE(cherk, CHERK)(const char *, const char *, int *, int *, const Kokkos::complex<float> *,
                                   const Kokkos::complex<float> *, int *, const Kokkos::complex<float> *,
                                   /* */ Kokkos::complex<float> *, int *);
void F77_BLAS_MANGLE(zherk, ZHERK)(const char *, const char *, int *, int *, const Kokkos::complex<double> *,
                                   const Kokkos::complex<double> *, int *, const Kokkos::complex<double> *,
                                   /* */ Kokkos::complex<double> *, int *);

///
/// Trsm
///

void F77_BLAS_MANGLE(strsm, STRSM)(const char *, const char *, const char *, const char *, int *, int *, const float *,
                                   const float *, int *,
                                   /* */ float *, int *);
void F77_BLAS_MANGLE(dtrsm, DTRSM)(const char *, const char *, const char *, const char *, int *, int *, const double *,
                                   const double *, int *,
                                   /* */ double *, int *);
void F77_BLAS_MANGLE(ctrsm, CTRSM)(const char *, const char *, const char *, const char *, int *, int *,
                                   const Kokkos::complex<float> *, const Kokkos::complex<float> *, int *,
                                   /* */ Kokkos::complex<float> *, int *);
void F77_BLAS_MANGLE(ztrsm, ZTRSM)(const char *, const char *, const char *, const char *, int *, int *,
                                   const Kokkos::complex<double> *, const Kokkos::complex<double> *, int *,
                                   /* */ Kokkos::complex<double> *, int *);
}

#define F77_FUNC_SGEMV F77_BLAS_MANGLE(sgemv, SGEMV)
#define F77_FUNC_DGEMV F77_BLAS_MANGLE(dgemv, DGEMV)
#define F77_FUNC_CGEMV F77_BLAS_MANGLE(cgemv, CGEMV)
#define F77_FUNC_ZGEMV F77_BLAS_MANGLE(zgemv, ZGEMV)

#define F77_FUNC_STRSV F77_BLAS_MANGLE(strsv, STRSV)
#define F77_FUNC_DTRSV F77_BLAS_MANGLE(dtrsv, DTRSV)
#define F77_FUNC_CTRSV F77_BLAS_MANGLE(ctrsv, CTRSV)
#define F77_FUNC_ZTRSV F77_BLAS_MANGLE(ztrsv, ZTRSV)

#define F77_FUNC_SGEMM F77_BLAS_MANGLE(sgemm, SGEMM)
#define F77_FUNC_DGEMM F77_BLAS_MANGLE(dgemm, DGEMM)
#define F77_FUNC_CGEMM F77_BLAS_MANGLE(cgemm, CGEMM)
#define F77_FUNC_ZGEMM F77_BLAS_MANGLE(zgemm, ZGEMM)

#define F77_FUNC_SSYRK F77_BLAS_MANGLE(ssyrk, SSYRK)
#define F77_FUNC_DSYRK F77_BLAS_MANGLE(dsyrk, DSYRK)
#define F77_FUNC_CHERK F77_BLAS_MANGLE(cherk, CHERK)
#define F77_FUNC_ZHERK F77_BLAS_MANGLE(zherk, ZHERK)

#define F77_FUNC_STRSM F77_BLAS_MANGLE(strsm, STRSM)
#define F77_FUNC_DTRSM F77_BLAS_MANGLE(dtrsm, DTRSM)
#define F77_FUNC_CTRSM F77_BLAS_MANGLE(ctrsm, CTRSM)
#define F77_FUNC_ZTRSM F77_BLAS_MANGLE(ztrsm, ZTRSM)

namespace Tacho {

///
/// float
///

template <>
int Blas<float>::gemv(const char trans, int m, int n, const float alpha, const float *a, int lda, const float *b,
                      int ldb, const float beta,
                      /* */ float *c, int ldc) {
  F77_FUNC_SGEMV(&trans, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<float>::gemv(cublasHandle_t handle, const cublasOperation_t trans, int m, int n, const float alpha,
                      const float *a, int lda, const float *b, int ldb, const float beta,
                      /* */ float *c, int ldc) {
  const int r_val = cublasSgemv(handle, trans, m, n, &alpha, a, lda, b, ldb, &beta, c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<float>::gemv(rocblas_handle handle, const rocblas_operation trans, int m, int n, const float alpha,
                      const float *a, int lda, const float *b, int ldb, const float beta,
                      /* */ float *c, int ldc) {
  const int r_val = rocblas_sgemv(handle, trans, m, n, &alpha, a, lda, b, ldb, &beta, c, ldc);
  return r_val;
}
#endif

template <>
int Blas<float>::trsv(const char uplo, const char transa, const char diag, int m, const float *a, int lda,
                      /* */ float *b, int ldb) {
  F77_FUNC_STRSV(&uplo, &transa, &diag, &m, a, &lda, b, &ldb);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<float>::trsv(cublasHandle_t handle, const cublasFillMode_t uplo, const cublasOperation_t transa,
                      const cublasDiagType_t diag, int m, const float *a, int lda,
                      /* */ float *b, int ldb) {
  const int r_val = cublasStrsv(handle, uplo, transa, diag, m, a, lda, b, ldb);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<float>::trsv(rocblas_handle handle, const rocblas_fill uplo, const rocblas_operation transa,
                      const rocblas_diagonal diag, int m, const float *a, int lda,
                      /* */ float *b, int ldb) {
  const int r_val = rocblas_strsv(handle, uplo, transa, diag, m, a, lda, b, ldb);
  return r_val;
}
#endif

template <>
int Blas<float>::gemm(const char transa, const char transb, int m, int n, int k, const float alpha, const float *a,
                      int lda, const float *b, int ldb, const float beta,
                      /* */ float *c, int ldc) {
  F77_FUNC_SGEMM(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<float>::gemm(cublasHandle_t handle, const cublasOperation_t transa, const cublasOperation_t transb, int m,
                      int n, int k, const float alpha, const float *a, int lda, const float *b, int ldb,
                      const float beta,
                      /* */ float *c, int ldc) {
  const int r_val = cublasSgemm(handle, transa, transb, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<float>::gemm(rocblas_handle handle, const rocblas_operation transa, const rocblas_operation transb, int m,
                      int n, int k, const float alpha, const float *a, int lda, const float *b, int ldb,
                      const float beta,
                      /* */ float *c, int ldc) {
  const int r_val = rocblas_sgemm(handle, transa, transb, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
  return r_val;
}
#endif

template <>
int Blas<float>::herk(const char uplo, const char trans, int n, int k, const float alpha, const float *a, int lda,
                      const float beta,
                      /* */ float *c, int ldc) {
  F77_FUNC_SSYRK(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<float>::herk(cublasHandle_t handle, const cublasFillMode_t uplo, const cublasOperation_t trans, int n, int k,
                      const float alpha, const float *a, int lda, const float beta,
                      /* */ float *c, int ldc) {
  const int r_val = cublasSsyrk(handle, uplo, trans, n, k, &alpha, a, lda, &beta, c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<float>::herk(rocblas_handle handle, const rocblas_fill uplo, const rocblas_operation trans, int n, int k,
                      const float alpha, const float *a, int lda, const float beta,
                      /* */ float *c, int ldc) {
  const int r_val = rocblas_ssyrk(handle, uplo, trans, n, k, &alpha, a, lda, &beta, c, ldc);
  return r_val;
}
#endif

template <>
int Blas<float>::trsm(const char side, const char uplo, const char transa, const char diag, int m, int n,
                      const float alpha, const float *a, int lda,
                      /* */ float *b, int ldb) {
  F77_FUNC_STRSM(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<float>::trsm(cublasHandle_t handle, const cublasSideMode_t side, const cublasFillMode_t uplo,
                      const cublasOperation_t transa, const cublasDiagType_t diag, int m, int n, const float alpha,
                      const float *a, int lda,
                      /* */ float *b, int ldb) {
  const int r_val = cublasStrsm(handle, side, uplo, transa, diag, m, n, &alpha, a, lda, b, ldb);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<float>::trsm(rocblas_handle handle, const rocblas_side side, const rocblas_fill uplo,
                      const rocblas_operation transa, const rocblas_diagonal diag, int m, int n, const float alpha,
                      const float *a, int lda,
                      /* */ float *b, int ldb) {
  const int r_val = rocblas_strsm(handle, side, uplo, transa, diag, m, n, &alpha, a, lda, b, ldb);
  return r_val;
}
#endif

///
/// double
///

template <>
int Blas<double>::gemv(const char trans, int m, int n, const double alpha, const double *a, int lda, const double *b,
                       int ldb, const double beta,
                       /* */ double *c, int ldc) {
  F77_FUNC_DGEMV(&trans, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<double>::gemv(cublasHandle_t handle, const cublasOperation_t trans, int m, int n, const double alpha,
                       const double *a, int lda, const double *b, int ldb, const double beta,
                       /* */ double *c, int ldc) {
  const int r_val = cublasDgemv(handle, trans, m, n, &alpha, a, lda, b, ldb, &beta, c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<double>::gemv(rocblas_handle handle, const rocblas_operation trans, int m, int n, const double alpha,
                       const double *a, int lda, const double *b, int ldb, const double beta,
                       /* */ double *c, int ldc) {
  const int r_val = rocblas_dgemv(handle, trans, m, n, &alpha, a, lda, b, ldb, &beta, c, ldc);
  return r_val;
}
#endif

template <>
int Blas<double>::trsv(const char uplo, const char transa, const char diag, int m, const double *a, int lda,
                       /* */ double *b, int ldb) {
  F77_FUNC_DTRSV(&uplo, &transa, &diag, &m, a, &lda, b, &ldb);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<double>::trsv(cublasHandle_t handle, const cublasFillMode_t uplo, const cublasOperation_t transa,
                       const cublasDiagType_t diag, int m, const double *a, int lda,
                       /* */ double *b, int ldb) {
  const int r_val = cublasDtrsv(handle, uplo, transa, diag, m, a, lda, b, ldb);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<double>::trsv(rocblas_handle handle, const rocblas_fill uplo, const rocblas_operation transa,
                       const rocblas_diagonal diag, int m, const double *a, int lda,
                       /* */ double *b, int ldb) {
  const int r_val = rocblas_dtrsv(handle, uplo, transa, diag, m, a, lda, b, ldb);
  return r_val;
}
#endif

template <>
int Blas<double>::gemm(const char transa, const char transb, int m, int n, int k, const double alpha, const double *a,
                       int lda, const double *b, int ldb, const double beta,
                       /* */ double *c, int ldc) {
  F77_FUNC_DGEMM(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<double>::gemm(cublasHandle_t handle, const cublasOperation_t transa, const cublasOperation_t transb, int m,
                       int n, int k, const double alpha, const double *a, int lda, const double *b, int ldb,
                       const double beta,
                       /* */ double *c, int ldc) {
  const int r_val = cublasDgemm(handle, transa, transb, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<double>::gemm(rocblas_handle handle, const rocblas_operation transa, const rocblas_operation transb, int m,
                       int n, int k, const double alpha, const double *a, int lda, const double *b, int ldb,
                       const double beta,
                       /* */ double *c, int ldc) {
  const int r_val = rocblas_dgemm(handle, transa, transb, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
  return r_val;
}
#endif

template <>
int Blas<double>::herk(const char uplo, const char trans, int n, int k, const double alpha, const double *a, int lda,
                       const double beta,
                       /* */ double *c, int ldc) {
  F77_FUNC_DSYRK(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<double>::herk(cublasHandle_t handle, const cublasFillMode_t uplo, const cublasOperation_t trans, int n, int k,
                       const double alpha, const double *a, int lda, const double beta,
                       /* */ double *c, int ldc) {
  const int r_val = cublasDsyrk(handle, uplo, trans, n, k, &alpha, a, lda, &beta, c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<double>::herk(rocblas_handle handle, const rocblas_fill uplo, const rocblas_operation trans, int n, int k,
                       const double alpha, const double *a, int lda, const double beta,
                       /* */ double *c, int ldc) {
  const int r_val = rocblas_dsyrk(handle, uplo, trans, n, k, &alpha, a, lda, &beta, c, ldc);
  return r_val;
}
#endif

template <>
int Blas<double>::trsm(const char side, const char uplo, const char transa, const char diag, int m, int n,
                       const double alpha, const double *a, int lda,
                       /* */ double *b, int ldb) {
  F77_FUNC_DTRSM(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<double>::trsm(cublasHandle_t handle, const cublasSideMode_t side, const cublasFillMode_t uplo,
                       const cublasOperation_t transa, const cublasDiagType_t diag, int m, int n, const double alpha,
                       const double *a, int lda,
                       /* */ double *b, int ldb) {
  const int r_val = cublasDtrsm(handle, side, uplo, transa, diag, m, n, &alpha, a, lda, b, ldb);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<double>::trsm(rocblas_handle handle, const rocblas_side side, const rocblas_fill uplo,
                       const rocblas_operation transa, const rocblas_diagonal diag, int m, int n, const double alpha,
                       const double *a, int lda,
                       /* */ double *b, int ldb) {
  const int r_val = rocblas_dtrsm(handle, side, uplo, transa, diag, m, n, &alpha, a, lda, b, ldb);
  return r_val;
}
#endif

///
/// Kokkos::complex<float>
///

template <>
int Blas<Kokkos::complex<float>>::gemv(const char trans, int m, int n, const Kokkos::complex<float> alpha,
                                       const Kokkos::complex<float> *a, int lda, const Kokkos::complex<float> *b,
                                       int ldb, const Kokkos::complex<float> beta,
                                       /* */ Kokkos::complex<float> *c, int ldc) {
  F77_FUNC_CGEMV(&trans, &m, &n, &alpha, (const Kokkos::complex<float> *)a, &lda, (const Kokkos::complex<float> *)b,
                 &ldb, &beta, (Kokkos::complex<float> *)c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<Kokkos::complex<float>>::gemv(cublasHandle_t handle, const cublasOperation_t trans, int m, int n,
                                       const Kokkos::complex<float> alpha, const Kokkos::complex<float> *a, int lda,
                                       const Kokkos::complex<float> *b, int ldb, const Kokkos::complex<float> beta,
                                       /* */ Kokkos::complex<float> *c, int ldc) {
  const int r_val = cublasCgemv(handle, trans, m, n, (const cuComplex *)&alpha, (const cuComplex *)a, lda,
                                (const cuComplex *)b, ldb, (const cuComplex *)&beta, (cuComplex *)c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<Kokkos::complex<float>>::gemv(rocblas_handle handle, const rocblas_operation trans, int m, int n,
                                       const Kokkos::complex<float> alpha, const Kokkos::complex<float> *a, int lda,
                                       const Kokkos::complex<float> *b, int ldb, const Kokkos::complex<float> beta,
                                       /* */ Kokkos::complex<float> *c, int ldc) {
  const int r_val = rocblas_cgemv(handle, trans, m, n, (const rocblas_float_complex *)&alpha,
                                  (const rocblas_float_complex *)a, lda, (const rocblas_float_complex *)b, ldb,
                                  (const rocblas_float_complex *)&beta, (rocblas_float_complex *)c, ldc);
  return r_val;
}
#endif

template <>
int Blas<Kokkos::complex<float>>::trsv(const char uplo, const char transa, const char diag, int m,
                                       const Kokkos::complex<float> *a, int lda,
                                       /* */ Kokkos::complex<float> *b, int ldb) {
  F77_FUNC_CTRSV(&uplo, &transa, &diag, &m, (const Kokkos::complex<float> *)a, &lda, (Kokkos::complex<float> *)b, &ldb);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<Kokkos::complex<float>>::trsv(cublasHandle_t handle, const cublasFillMode_t uplo,
                                       const cublasOperation_t transa, const cublasDiagType_t diag, int m,
                                       const Kokkos::complex<float> *a, int lda,
                                       /* */ Kokkos::complex<float> *b, int ldb) {
  const int r_val = cublasCtrsv(handle, uplo, transa, diag, m, (const cuComplex *)a, lda, (cuComplex *)b, ldb);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<Kokkos::complex<float>>::trsv(rocblas_handle handle, const rocblas_fill uplo, const rocblas_operation transa,
                                       const rocblas_diagonal diag, int m, const Kokkos::complex<float> *a, int lda,
                                       /* */ Kokkos::complex<float> *b, int ldb) {
  const int r_val = rocblas_ctrsv(handle, uplo, transa, diag, m, (const rocblas_float_complex *)a, lda,
                                  (rocblas_float_complex *)b, ldb);
  return r_val;
}
#endif

template <>
int Blas<Kokkos::complex<float>>::gemm(const char transa, const char transb, int m, int n, int k,
                                       const Kokkos::complex<float> alpha, const Kokkos::complex<float> *a, int lda,
                                       const Kokkos::complex<float> *b, int ldb, const Kokkos::complex<float> beta,
                                       /* */ Kokkos::complex<float> *c, int ldc) {
  F77_FUNC_CGEMM(&transa, &transb, &m, &n, &k, &alpha, (const Kokkos::complex<float> *)a, &lda,
                 (const Kokkos::complex<float> *)b, &ldb, &beta, (Kokkos::complex<float> *)c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<Kokkos::complex<float>>::gemm(cublasHandle_t handle, const cublasOperation_t transa,
                                       const cublasOperation_t transb, int m, int n, int k,
                                       const Kokkos::complex<float> alpha, const Kokkos::complex<float> *a, int lda,
                                       const Kokkos::complex<float> *b, int ldb, const Kokkos::complex<float> beta,
                                       /* */ Kokkos::complex<float> *c, int ldc) {
  const int r_val = cublasCgemm(handle, transa, transb, m, n, k, (const cuComplex *)&alpha, (const cuComplex *)a, lda,
                                (const cuComplex *)b, ldb, (const cuComplex *)&beta, (cuComplex *)c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<Kokkos::complex<float>>::gemm(rocblas_handle handle, const rocblas_operation transa,
                                       const rocblas_operation transb, int m, int n, int k,
                                       const Kokkos::complex<float> alpha, const Kokkos::complex<float> *a, int lda,
                                       const Kokkos::complex<float> *b, int ldb, const Kokkos::complex<float> beta,
                                       /* */ Kokkos::complex<float> *c, int ldc) {
  const int r_val = rocblas_cgemm(handle, transa, transb, m, n, k, (const rocblas_float_complex *)&alpha,
                                  (const rocblas_float_complex *)a, lda, (const rocblas_float_complex *)b, ldb,
                                  (const rocblas_float_complex *)&beta, (rocblas_float_complex *)c, ldc);
  return r_val;
}
#endif

template <>
int Blas<Kokkos::complex<float>>::herk(const char uplo, const char trans, int n, int k,
                                       const Kokkos::complex<float> alpha, const Kokkos::complex<float> *a, int lda,
                                       const Kokkos::complex<float> beta,
                                       /* */ Kokkos::complex<float> *c, int ldc) {
  F77_FUNC_CHERK(&uplo, &trans, &n, &k, &alpha, (const Kokkos::complex<float> *)a, &lda, &beta,
                 (Kokkos::complex<float> *)c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<Kokkos::complex<float>>::herk(cublasHandle_t handle, const cublasFillMode_t uplo,
                                       const cublasOperation_t trans, int n, int k, const Kokkos::complex<float> alpha,
                                       const Kokkos::complex<float> *a, int lda, const Kokkos::complex<float> beta,
                                       /* */ Kokkos::complex<float> *c, int ldc) {
  const int r_val = cublasCherk(handle, uplo, trans, n, k, (const float *)&alpha, (const cuComplex *)a, lda,
                                (const float *)&beta, (cuComplex *)c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<Kokkos::complex<float>>::herk(rocblas_handle handle, const rocblas_fill uplo, const rocblas_operation trans,
                                       int n, int k, const Kokkos::complex<float> alpha,
                                       const Kokkos::complex<float> *a, int lda, const Kokkos::complex<float> beta,
                                       /* */ Kokkos::complex<float> *c, int ldc) {
  const int r_val = rocblas_cherk(handle, uplo, trans, n, k, (const float *)&alpha, (const rocblas_float_complex *)a,
                                  lda, (const float *)&beta, (rocblas_float_complex *)c, ldc);
  return r_val;
}
#endif

template <>
int Blas<Kokkos::complex<float>>::trsm(const char side, const char uplo, const char transa, const char diag, int m,
                                       int n, const Kokkos::complex<float> alpha, const Kokkos::complex<float> *a,
                                       int lda,
                                       /* */ Kokkos::complex<float> *b, int ldb) {
  F77_FUNC_CTRSM(&side, &uplo, &transa, &diag, &m, &n, &alpha, (const Kokkos::complex<float> *)a, &lda,
                 (Kokkos::complex<float> *)b, &ldb);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<Kokkos::complex<float>>::trsm(cublasHandle_t handle, const cublasSideMode_t side, const cublasFillMode_t uplo,
                                       const cublasOperation_t transa, const cublasDiagType_t diag, int m, int n,
                                       const Kokkos::complex<float> alpha, const Kokkos::complex<float> *a, int lda,
                                       /* */ Kokkos::complex<float> *b, int ldb) {
  const int r_val = cublasCtrsm(handle, side, uplo, transa, diag, m, n, (const cuComplex *)&alpha, (const cuComplex *)a,
                                lda, (cuComplex *)b, ldb);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<Kokkos::complex<float>>::trsm(rocblas_handle handle, const rocblas_side side, const rocblas_fill uplo,
                                       const rocblas_operation transa, const rocblas_diagonal diag, int m, int n,
                                       const Kokkos::complex<float> alpha, const Kokkos::complex<float> *a, int lda,
                                       /* */ Kokkos::complex<float> *b, int ldb) {
  const int r_val = rocblas_ctrsm(handle, side, uplo, transa, diag, m, n, (const rocblas_float_complex *)&alpha,
                                  (const rocblas_float_complex *)a, lda, (rocblas_float_complex *)b, ldb);
  return r_val;
}
#endif

///
/// Kokkos::complex<double>
///

template <>
int Blas<Kokkos::complex<double>>::gemv(const char trans, int m, int n, const Kokkos::complex<double> alpha,
                                        const Kokkos::complex<double> *a, int lda, const Kokkos::complex<double> *b,
                                        int ldb, const Kokkos::complex<double> beta,
                                        /* */ Kokkos::complex<double> *c, int ldc) {
  F77_FUNC_ZGEMV(&trans, &m, &n, &alpha, (const Kokkos::complex<double> *)a, &lda, (const Kokkos::complex<double> *)b,
                 &ldb, &beta, (Kokkos::complex<double> *)c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<Kokkos::complex<double>>::gemv(cublasHandle_t handle, const cublasOperation_t trans, int m, int n,
                                        const Kokkos::complex<double> alpha, const Kokkos::complex<double> *a, int lda,
                                        const Kokkos::complex<double> *b, int ldb, const Kokkos::complex<double> beta,
                                        /* */ Kokkos::complex<double> *c, int ldc) {
  const int r_val =
      cublasZgemv(handle, trans, m, n, (const cuDoubleComplex *)&alpha, (const cuDoubleComplex *)a, lda,
                  (const cuDoubleComplex *)b, ldb, (const cuDoubleComplex *)&beta, (cuDoubleComplex *)c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<Kokkos::complex<double>>::gemv(rocblas_handle handle, const rocblas_operation trans, int m, int n,
                                        const Kokkos::complex<double> alpha, const Kokkos::complex<double> *a, int lda,
                                        const Kokkos::complex<double> *b, int ldb, const Kokkos::complex<double> beta,
                                        /* */ Kokkos::complex<double> *c, int ldc) {
  const int r_val = rocblas_zgemv(handle, trans, m, n, (const rocblas_double_complex *)&alpha,
                                  (const rocblas_double_complex *)a, lda, (const rocblas_double_complex *)b, ldb,
                                  (const rocblas_double_complex *)&beta, (rocblas_double_complex *)c, ldc);
  return r_val;
}
#endif

template <>
int Blas<Kokkos::complex<double>>::trsv(const char uplo, const char transa, const char diag, int m,
                                        const Kokkos::complex<double> *a, int lda,
                                        /* */ Kokkos::complex<double> *b, int ldb) {
  F77_FUNC_ZTRSV(&uplo, &transa, &diag, &m, (const Kokkos::complex<double> *)a, &lda, (Kokkos::complex<double> *)b,
                 &ldb);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<Kokkos::complex<double>>::trsv(cublasHandle_t handle, const cublasFillMode_t uplo,
                                        const cublasOperation_t transa, const cublasDiagType_t diag, int m,
                                        const Kokkos::complex<double> *a, int lda,
                                        /* */ Kokkos::complex<double> *b, int ldb) {
  const int r_val =
      cublasZtrsv(handle, uplo, transa, diag, m, (const cuDoubleComplex *)a, lda, (cuDoubleComplex *)b, ldb);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<Kokkos::complex<double>>::trsv(rocblas_handle handle, const rocblas_fill uplo, const rocblas_operation transa,
                                        const rocblas_diagonal diag, int m, const Kokkos::complex<double> *a, int lda,
                                        /* */ Kokkos::complex<double> *b, int ldb) {
  const int r_val = rocblas_ztrsv(handle, uplo, transa, diag, m, (const rocblas_double_complex *)a, lda,
                                  (rocblas_double_complex *)b, ldb);
  return r_val;
}
#endif

template <>
int Blas<Kokkos::complex<double>>::gemm(const char transa, const char transb, int m, int n, int k,
                                        const Kokkos::complex<double> alpha, const Kokkos::complex<double> *a, int lda,
                                        const Kokkos::complex<double> *b, int ldb, const Kokkos::complex<double> beta,
                                        /* */ Kokkos::complex<double> *c, int ldc) {
  F77_FUNC_ZGEMM(&transa, &transb, &m, &n, &k, &alpha, (const Kokkos::complex<double> *)a, &lda,
                 (const Kokkos::complex<double> *)b, &ldb, &beta, (Kokkos::complex<double> *)c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<Kokkos::complex<double>>::gemm(cublasHandle_t handle, const cublasOperation_t transa,
                                        const cublasOperation_t transb, int m, int n, int k,
                                        const Kokkos::complex<double> alpha, const Kokkos::complex<double> *a, int lda,
                                        const Kokkos::complex<double> *b, int ldb, const Kokkos::complex<double> beta,
                                        /* */ Kokkos::complex<double> *c, int ldc) {
  const int r_val =
      cublasZgemm(handle, transa, transb, m, n, k, (const cuDoubleComplex *)&alpha, (const cuDoubleComplex *)a, lda,
                  (const cuDoubleComplex *)b, ldb, (const cuDoubleComplex *)&beta, (cuDoubleComplex *)c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<Kokkos::complex<double>>::gemm(rocblas_handle handle, const rocblas_operation transa,
                                        const rocblas_operation transb, int m, int n, int k,
                                        const Kokkos::complex<double> alpha, const Kokkos::complex<double> *a, int lda,
                                        const Kokkos::complex<double> *b, int ldb, const Kokkos::complex<double> beta,
                                        /* */ Kokkos::complex<double> *c, int ldc) {
  const int r_val = rocblas_zgemm(handle, transa, transb, m, n, k, (const rocblas_double_complex *)&alpha,
                                  (const rocblas_double_complex *)a, lda, (const rocblas_double_complex *)b, ldb,
                                  (const rocblas_double_complex *)&beta, (rocblas_double_complex *)c, ldc);
  return r_val;
}
#endif

template <>
int Blas<Kokkos::complex<double>>::herk(const char uplo, const char trans, int n, int k,
                                        const Kokkos::complex<double> alpha, const Kokkos::complex<double> *a, int lda,
                                        const Kokkos::complex<double> beta,
                                        /* */ Kokkos::complex<double> *c, int ldc) {
  F77_FUNC_ZHERK(&uplo, &trans, &n, &k, &alpha, (const Kokkos::complex<double> *)a, &lda, &beta,
                 (Kokkos::complex<double> *)c, &ldc);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<Kokkos::complex<double>>::herk(cublasHandle_t handle, const cublasFillMode_t uplo,
                                        const cublasOperation_t trans, int n, int k,
                                        const Kokkos::complex<double> alpha, const Kokkos::complex<double> *a, int lda,
                                        const Kokkos::complex<double> beta,
                                        /* */ Kokkos::complex<double> *c, int ldc) {
  const int r_val = cublasZherk(handle, uplo, trans, n, k, (const double *)&alpha, (const cuDoubleComplex *)a, lda,
                                (const double *)&beta, (cuDoubleComplex *)c, ldc);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<Kokkos::complex<double>>::herk(rocblas_handle handle, const rocblas_fill uplo, const rocblas_operation trans,
                                        int n, int k, const Kokkos::complex<double> alpha,
                                        const Kokkos::complex<double> *a, int lda, const Kokkos::complex<double> beta,
                                        /* */ Kokkos::complex<double> *c, int ldc) {
  const int r_val = rocblas_zherk(handle, uplo, trans, n, k, (const double *)&alpha, (const rocblas_double_complex *)a,
                                  lda, (const double *)&beta, (rocblas_double_complex *)c, ldc);
  return r_val;
}
#endif

template <>
int Blas<Kokkos::complex<double>>::trsm(const char side, const char uplo, const char transa, const char diag, int m,
                                        int n, const Kokkos::complex<double> alpha, const Kokkos::complex<double> *a,
                                        int lda,
                                        /* */ Kokkos::complex<double> *b, int ldb) {
  F77_FUNC_ZTRSM(&side, &uplo, &transa, &diag, &m, &n, &alpha, (const Kokkos::complex<double> *)a, &lda,
                 (Kokkos::complex<double> *)b, &ldb);
  return 0;
}
#if defined(TACHO_ENABLE_CUBLAS)
template <>
int Blas<Kokkos::complex<double>>::trsm(cublasHandle_t handle, const cublasSideMode_t side, const cublasFillMode_t uplo,
                                        const cublasOperation_t transa, const cublasDiagType_t diag, int m, int n,
                                        const Kokkos::complex<double> alpha, const Kokkos::complex<double> *a, int lda,
                                        /* */ Kokkos::complex<double> *b, int ldb) {
  const int r_val = cublasZtrsm(handle, side, uplo, transa, diag, m, n, (const cuDoubleComplex *)&alpha,
                                (const cuDoubleComplex *)a, lda, (cuDoubleComplex *)b, ldb);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCBLAS)
template <>
int Blas<Kokkos::complex<double>>::trsm(rocblas_handle handle, const rocblas_side side, const rocblas_fill uplo,
                                        const rocblas_operation transa, const rocblas_diagonal diag, int m, int n,
                                        const Kokkos::complex<double> alpha, const Kokkos::complex<double> *a, int lda,
                                        /* */ Kokkos::complex<double> *b, int ldb) {
  const int r_val = rocblas_ztrsm(handle, side, uplo, transa, diag, m, n, (const rocblas_double_complex *)&alpha,
                                  (const rocblas_double_complex *)a, lda, (rocblas_double_complex *)b, ldb);
  return r_val;
}
#endif

///
/// std::complex<float>
///

template <>
int Blas<std::complex<float>>::gemv(const char trans, int m, int n, const std::complex<float> alpha,
                                    const std::complex<float> *a, int lda, const std::complex<float> *b, int ldb,
                                    const std::complex<float> beta,
                                    /* */ std::complex<float> *c, int ldc) {
  F77_FUNC_CGEMV(&trans, &m, &n, (const Kokkos::complex<float> *)&alpha, (const Kokkos::complex<float> *)a, &lda,
                 (const Kokkos::complex<float> *)b, &ldb, (const Kokkos::complex<float> *)&beta,
                 (Kokkos::complex<float> *)c, &ldc);
  return 0;
}
template <>
int Blas<std::complex<float>>::trsv(const char uplo, const char transa, const char diag, int m,
                                    const std::complex<float> *a, int lda,
                                    /* */ std::complex<float> *b, int ldb) {
  F77_FUNC_CTRSV(&uplo, &transa, &diag, &m, (const Kokkos::complex<float> *)a, &lda, (Kokkos::complex<float> *)b, &ldb);
  return 0;
}
template <>
int Blas<std::complex<float>>::gemm(const char transa, const char transb, int m, int n, int k,
                                    const std::complex<float> alpha, const std::complex<float> *a, int lda,
                                    const std::complex<float> *b, int ldb, const std::complex<float> beta,
                                    /* */ std::complex<float> *c, int ldc) {
  F77_FUNC_CGEMM(&transa, &transb, &m, &n, &k, (const Kokkos::complex<float> *)&alpha,
                 (const Kokkos::complex<float> *)a, &lda, (const Kokkos::complex<float> *)b, &ldb,
                 (const Kokkos::complex<float> *)&beta, (Kokkos::complex<float> *)c, &ldc);
  return 0;
}
template <>
int Blas<std::complex<float>>::herk(const char transa, const char transb, int n, int k, const std::complex<float> alpha,
                                    const std::complex<float> *a, int lda, const std::complex<float> beta,
                                    /* */ std::complex<float> *c, int ldc) {
  F77_FUNC_CHERK(&transa, &transb, &n, &k, (const Kokkos::complex<float> *)&alpha, (const Kokkos::complex<float> *)a,
                 &lda, (const Kokkos::complex<float> *)&beta, (Kokkos::complex<float> *)c, &ldc);
  return 0;
}
template <>
int Blas<std::complex<float>>::trsm(const char side, const char uplo, const char transa, const char diag, int m, int n,
                                    const std::complex<float> alpha, const std::complex<float> *a, int lda,
                                    /* */ std::complex<float> *b, int ldb) {
  F77_FUNC_CTRSM(&side, &uplo, &transa, &diag, &m, &n, (const Kokkos::complex<float> *)&alpha,
                 (const Kokkos::complex<float> *)a, &lda, (Kokkos::complex<float> *)b, &ldb);
  return 0;
}

///
/// std::complex<double>
///

template <>
int Blas<std::complex<double>>::gemv(const char trans, int m, int n, const std::complex<double> alpha,
                                     const std::complex<double> *a, int lda, const std::complex<double> *b, int ldb,
                                     const std::complex<double> beta,
                                     /* */ std::complex<double> *c, int ldc) {
  F77_FUNC_ZGEMV(&trans, &m, &n, (const Kokkos::complex<double> *)&alpha, (const Kokkos::complex<double> *)a, &lda,
                 (const Kokkos::complex<double> *)b, &ldb, (const Kokkos::complex<double> *)&beta,
                 (Kokkos::complex<double> *)c, &ldc);
  return 0;
}
template <>
int Blas<std::complex<double>>::trsv(const char uplo, const char transa, const char diag, int m,
                                     const std::complex<double> *a, int lda,
                                     /* */ std::complex<double> *b, int ldb) {
  F77_FUNC_ZTRSV(&uplo, &transa, &diag, &m, (const Kokkos::complex<double> *)a, &lda, (Kokkos::complex<double> *)b,
                 &ldb);
  return 0;
}
template <>
int Blas<std::complex<double>>::gemm(const char transa, const char transb, int m, int n, int k,
                                     const std::complex<double> alpha, const std::complex<double> *a, int lda,
                                     const std::complex<double> *b, int ldb, const std::complex<double> beta,
                                     /* */ std::complex<double> *c, int ldc) {
  F77_FUNC_ZGEMM(&transa, &transb, &m, &n, &k, (const Kokkos::complex<double> *)&alpha,
                 (const Kokkos::complex<double> *)a, &lda, (const Kokkos::complex<double> *)b, &ldb,
                 (const Kokkos::complex<double> *)&beta, (Kokkos::complex<double> *)c, &ldc);
  return 0;
}
template <>
int Blas<std::complex<double>>::herk(const char transa, const char transb, int n, int k,
                                     const std::complex<double> alpha, const std::complex<double> *a, int lda,
                                     const std::complex<double> beta,
                                     /* */ std::complex<double> *c, int ldc) {
  F77_FUNC_ZHERK(&transa, &transb, &n, &k, (const Kokkos::complex<double> *)&alpha, (const Kokkos::complex<double> *)a,
                 &lda, (const Kokkos::complex<double> *)&beta, (Kokkos::complex<double> *)c, &ldc);
  return 0;
}
template <>
int Blas<std::complex<double>>::trsm(const char side, const char uplo, const char transa, const char diag, int m, int n,
                                     const std::complex<double> alpha, const std::complex<double> *a, int lda,
                                     /* */ std::complex<double> *b, int ldb) {
  F77_FUNC_ZTRSM(&side, &uplo, &transa, &diag, &m, &n, (const Kokkos::complex<double> *)&alpha,
                 (const Kokkos::complex<double> *)a, &lda, (Kokkos::complex<double> *)b, &ldb);
  return 0;
}

} // namespace Tacho
