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
/// \file  Tacho_Lapack_External.hpp
/// \brief Lapack wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Lapack_External.hpp"
#include "Tacho_config.h"

#include "Kokkos_Core.hpp"

extern "C" {
///
/// Cholesky
///
void F77_BLAS_MANGLE(spotrf, SPOTRF)(const char *, const int *, float *, const int *, int *);
void F77_BLAS_MANGLE(dpotrf, DPOTRF)(const char *, const int *, double *, const int *, int *);
void F77_BLAS_MANGLE(cpotrf, CPOTRF)(const char *, const int *, Kokkos::complex<float> *, const int *, int *);
void F77_BLAS_MANGLE(zpotrf, ZPOTRF)(const char *, const int *, Kokkos::complex<double> *, const int *, int *);

///
/// LDLt
///
void F77_BLAS_MANGLE(ssytrf, SSYTRF)(const char *, const int *, float *, const int *, int *, float *, int *, int *);
void F77_BLAS_MANGLE(dsytrf, DSYTRF)(const char *, const int *, double *, const int *, int *, double *, int *, int *);
void F77_BLAS_MANGLE(csytrf, CSYTRF)(const char *, const int *, Kokkos::complex<float> *, const int *, int *,
                                     Kokkos::complex<float> *, int *, int *);
void F77_BLAS_MANGLE(zsytrf, ZSYTRF)(const char *, const int *, Kokkos::complex<double> *, const int *, int *,
                                     Kokkos::complex<double> *, int *, int *);

///
/// LU
/// M, N, A, LDA, IPIV, INFO )
void F77_BLAS_MANGLE(sgetrf, SGETRF)(const int *, const int *, float *, const int *, int *, int *);
void F77_BLAS_MANGLE(dgetrf, DGETRF)(const int *, const int *, double *, const int *, int *, int *);
void F77_BLAS_MANGLE(cgetrf, CGETRF)(const int *, const int *, Kokkos::complex<float> *, const int *, int *, int *);
void F77_BLAS_MANGLE(zgetrf, ZGETRF)(const int *, const int *, Kokkos::complex<double> *, const int *, int *, int *);
}

namespace Tacho {

#define F77_FUNC_SPOTRF F77_BLAS_MANGLE(spotrf, SPOTRF)
#define F77_FUNC_DPOTRF F77_BLAS_MANGLE(dpotrf, DPOTRF)
#define F77_FUNC_CPOTRF F77_BLAS_MANGLE(cpotrf, CPOTRF)
#define F77_FUNC_ZPOTRF F77_BLAS_MANGLE(zpotrf, ZPOTRF)

#define F77_FUNC_SSYTRF F77_BLAS_MANGLE(ssytrf, SSYTRF)
#define F77_FUNC_DSYTRF F77_BLAS_MANGLE(dsytrf, DSYTRF)
#define F77_FUNC_CSYTRF F77_BLAS_MANGLE(csytrf, CSYTRF)
#define F77_FUNC_ZSYTRF F77_BLAS_MANGLE(zsytrf, ZSYTRF)

#define F77_FUNC_SGETRF F77_BLAS_MANGLE(sgetrf, SGETRF)
#define F77_FUNC_DGETRF F77_BLAS_MANGLE(dgetrf, DGETRF)
#define F77_FUNC_CGETRF F77_BLAS_MANGLE(cgetrf, CGETRF)
#define F77_FUNC_ZGETRF F77_BLAS_MANGLE(zgetrf, ZGETRF)

template <> int Lapack<float>::potrf(const char uplo, const int m, float *a, const int lda, int *info) {
  F77_FUNC_SPOTRF(&uplo, &m, a, &lda, info);
  return 0;
}
#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<float>::potrf_buffersize(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m, float *a,
                                    const int lda, int *lwork) {
  const int r_val = cusolverDnSpotrf_bufferSize(handle, uplo, m, a, lda, lwork);
  return r_val;
}

template <>
int Lapack<float>::potrf(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m, float *a, const int lda,
                         float *w, const int lwork, int *dev) {
  const int r_val = cusolverDnSpotrf(handle, uplo, m, a, lda, w, lwork, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<float>::potrf(rocblas_handle handle, const rocblas_fill uplo, const int m, float *a, const int lda,
                         int *dev) {
  const int r_val = rocsolver_spotrf(handle, uplo, m, a, lda, dev);
  return r_val;
}
#endif

template <>
int Lapack<float>::sytrf(const char uplo, const int m, float *a, const int lda, int *ipiv, float *work, int lwork,
                         int *info) {
  F77_FUNC_SSYTRF(&uplo, &m, a, &lda, ipiv, work, &lwork, info);
  return 0;
}

#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<float>::sytrf_buffersize(cusolverDnHandle_t handle, const int m, float *a, const int lda, int *lwork) {
  const int r_val = cusolverDnSsytrf_bufferSize(handle, m, a, lda, lwork);
  return r_val;
}

template <>
int Lapack<float>::sytrf(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m, float *a, const int lda,
                         int *ipiv, float *w, const int lwork, int *dev) {
  const int r_val = cusolverDnSsytrf(handle, uplo, m, a, lda, ipiv, w, lwork, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<float>::sytrf(rocblas_handle handle, const rocblas_fill uplo, const int m, float *a, const int lda,
                         int *ipiv, int *dev) {
#if defined(TACHO_ENABLE_ROCSOLVER_SYTRF)
  const int r_val = rocsolver_ssytrf(handle, uplo, m, a, lda, ipiv, dev);
  return r_val;
#else
  throw std::logic_error("Error: sytrf is not avail in old ROCM");
  return -1;
#endif
}
#endif

template <> int Lapack<float>::getrf(const int m, const int n, float *a, const int lda, int *ipiv, int *info) {
  F77_FUNC_SGETRF(&m, &n, a, &lda, ipiv, info);
  return 0;
}

#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<float>::getrf_buffersize(cusolverDnHandle_t handle, const int m, const int n, float *a, const int lda,
                                    int *lwork) {
  const int r_val = cusolverDnSgetrf_bufferSize(handle, m, n, a, lda, lwork);
  return r_val;
}

template <>
int Lapack<float>::getrf(cusolverDnHandle_t handle, const int m, const int n, float *a, const int lda, float *w,
                         int *ipiv, int *dev) {
  const int r_val = cusolverDnSgetrf(handle, m, n, a, lda, w, ipiv, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<float>::getrf(rocblas_handle handle, const int m, const int n, float *a, const int lda, int *ipiv,
                         int *dev) {
  const int r_val = rocsolver_sgetrf(handle, m, n, a, lda, ipiv, dev);
  return r_val;
}
#endif

template <> int Lapack<double>::potrf(const char uplo, const int m, double *a, const int lda, int *info) {
  F77_FUNC_DPOTRF(&uplo, &m, a, &lda, info);
  return 0;
}
#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<double>::potrf_buffersize(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m, double *a,
                                     const int lda, int *lwork) {
  const int r_val = cusolverDnDpotrf_bufferSize(handle, uplo, m, a, lda, lwork);
  return r_val;
}

template <>
int Lapack<double>::potrf(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m, double *a, const int lda,
                          double *w, const int lwork, int *dev) {
  const int r_val = cusolverDnDpotrf(handle, uplo, m, a, lda, w, lwork, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<double>::potrf(rocblas_handle handle, const rocblas_fill uplo, const int m, double *a, const int lda,
                          int *dev) {
  const int r_val = rocsolver_dpotrf(handle, uplo, m, a, lda, dev);
  return r_val;
}
#endif

template <>
int Lapack<double>::sytrf(const char uplo, const int m, double *a, const int lda, int *ipiv, double *work, int lwork,
                          int *info) {
  F77_FUNC_DSYTRF(&uplo, &m, a, &lda, ipiv, work, &lwork, info);
  return 0;
}
#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<double>::sytrf_buffersize(cusolverDnHandle_t handle, const int m, double *a, const int lda, int *lwork) {
  const int r_val = cusolverDnDsytrf_bufferSize(handle, m, a, lda, lwork);
  return r_val;
}

template <>
int Lapack<double>::sytrf(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m, double *a, const int lda,
                          int *ipiv, double *w, const int lwork, int *dev) {
  const int r_val = cusolverDnDsytrf(handle, uplo, m, a, lda, ipiv, w, lwork, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<double>::sytrf(rocblas_handle handle, const rocblas_fill uplo, const int m, double *a, const int lda,
                          int *ipiv, int *dev) {
#if defined(TACHO_ENABLE_ROCSOLVER_SYTRF)
  const int r_val = rocsolver_dsytrf(handle, uplo, m, a, lda, ipiv, dev);
  return r_val;
#else
  throw std::logic_error("Error: sytrf is not avail in old ROCM");
  return -1;
#endif
}
#endif

template <> int Lapack<double>::getrf(const int m, const int n, double *a, const int lda, int *ipiv, int *info) {
  F77_FUNC_DGETRF(&m, &n, a, &lda, ipiv, info);
  return 0;
}

#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<double>::getrf_buffersize(cusolverDnHandle_t handle, const int m, const int n, double *a, const int lda,
                                     int *lwork) {
  const int r_val = cusolverDnDgetrf_bufferSize(handle, m, n, a, lda, lwork);
  return r_val;
}

template <>
int Lapack<double>::getrf(cusolverDnHandle_t handle, const int m, const int n, double *a, const int lda, double *w,
                          int *ipiv, int *dev) {
  const int r_val = cusolverDnDgetrf(handle, m, n, a, lda, w, ipiv, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<double>::getrf(rocblas_handle handle, const int m, const int n, double *a, const int lda, int *ipiv,
                          int *dev) {
  const int r_val = rocsolver_dgetrf(handle, m, n, a, lda, ipiv, dev);
  return r_val;
}
#endif

template <>
int Lapack<Kokkos::complex<float>>::potrf(const char uplo, const int m, Kokkos::complex<float> *a, const int lda,
                                          int *info) {
  F77_FUNC_CPOTRF(&uplo, &m, a, &lda, info);
  return 0;
}
#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<Kokkos::complex<float>>::potrf_buffersize(cusolverDnHandle_t handle, const cublasFillMode_t uplo,
                                                     const int m, Kokkos::complex<float> *a, const int lda,
                                                     int *lwork) {
  const int r_val = cusolverDnCpotrf_bufferSize(handle, uplo, m, (cuComplex *)a, lda, lwork);
  return r_val;
}

template <>
int Lapack<Kokkos::complex<float>>::potrf(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m,
                                          Kokkos::complex<float> *a, const int lda, Kokkos::complex<float> *w,
                                          const int lwork, int *dev) {
  const int r_val = cusolverDnCpotrf(handle, uplo, m, (cuComplex *)a, lda, (cuComplex *)w, lwork, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<Kokkos::complex<float>>::potrf(rocblas_handle handle, const rocblas_fill uplo, const int m,
                                          Kokkos::complex<float> *a, const int lda, int *dev) {
  const int r_val = rocsolver_cpotrf(handle, uplo, m, (rocblas_float_complex *)a, lda, dev);
  return r_val;
}
#endif

template <>
int Lapack<Kokkos::complex<float>>::sytrf(const char uplo, const int m, Kokkos::complex<float> *a, const int lda,
                                          int *ipiv, Kokkos::complex<float> *work, int lwork, int *info) {
  F77_FUNC_CSYTRF(&uplo, &m, a, &lda, ipiv, work, &lwork, info);
  return 0;
}

#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<Kokkos::complex<float>>::sytrf_buffersize(cusolverDnHandle_t handle, const int m, Kokkos::complex<float> *a,
                                                     const int lda, int *lwork) {
  const int r_val = cusolverDnCsytrf_bufferSize(handle, m, (cuComplex *)a, lda, lwork);
  return r_val;
}

template <>
int Lapack<Kokkos::complex<float>>::sytrf(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m,
                                          Kokkos::complex<float> *a, const int lda, int *ipiv,
                                          Kokkos::complex<float> *w, const int lwork, int *dev) {
  const int r_val = cusolverDnCsytrf(handle, uplo, m, (cuComplex *)a, lda, ipiv, (cuComplex *)w, lwork, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<Kokkos::complex<float>>::sytrf(rocblas_handle handle, const rocblas_fill uplo, const int m,
                                          Kokkos::complex<float> *a, const int lda, int *ipiv, int *dev) {
#if defined(TACHO_ENABLE_ROCSOLVER_SYTRF)
  const int r_val = rocsolver_csytrf(handle, uplo, m, (rocblas_float_complex *)a, lda, ipiv, dev);
  return r_val;
#else
  throw std::logic_error("Error: sytrf is not avail in old ROCM");
  return -1;
#endif
}
#endif

template <>
int Lapack<Kokkos::complex<float>>::getrf(const int m, const int n, Kokkos::complex<float> *a, const int lda, int *ipiv,
                                          int *info) {
  F77_FUNC_CGETRF(&m, &n, a, &lda, ipiv, info);
  return 0;
}

#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<Kokkos::complex<float>>::getrf_buffersize(cusolverDnHandle_t handle, const int m, const int n,
                                                     Kokkos::complex<float> *a, const int lda, int *lwork) {
  const int r_val = cusolverDnCgetrf_bufferSize(handle, m, n, (cuComplex *)a, lda, lwork);
  return r_val;
}

template <>
int Lapack<Kokkos::complex<float>>::getrf(cusolverDnHandle_t handle, const int m, const int n,
                                          Kokkos::complex<float> *a, const int lda, Kokkos::complex<float> *w,
                                          int *ipiv, int *dev) {
  const int r_val = cusolverDnCgetrf(handle, m, n, (cuComplex *)a, lda, (cuComplex *)w, ipiv, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<Kokkos::complex<float>>::getrf(rocblas_handle handle, const int m, const int n, Kokkos::complex<float> *a,
                                          const int lda, int *ipiv, int *dev) {
  const int r_val = rocsolver_cgetrf(handle, m, n, (rocblas_float_complex *)a, lda, ipiv, dev);
  return r_val;
}
#endif

template <>
int Lapack<Kokkos::complex<double>>::potrf(const char uplo, const int m, Kokkos::complex<double> *a, const int lda,
                                           int *info) {
  F77_FUNC_ZPOTRF(&uplo, &m, a, &lda, info);
  return 0;
}
#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<Kokkos::complex<double>>::potrf_buffersize(cusolverDnHandle_t handle, const cublasFillMode_t uplo,
                                                      const int m, Kokkos::complex<double> *a, const int lda,
                                                      int *lwork) {
  const int r_val = cusolverDnZpotrf_bufferSize(handle, uplo, m, (cuDoubleComplex *)a, lda, lwork);
  return r_val;
}

template <>
int Lapack<Kokkos::complex<double>>::potrf(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m,
                                           Kokkos::complex<double> *a, const int lda, Kokkos::complex<double> *w,
                                           const int lwork, int *dev) {
  const int r_val = cusolverDnZpotrf(handle, uplo, m, (cuDoubleComplex *)a, lda, (cuDoubleComplex *)w, lwork, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<Kokkos::complex<double>>::potrf(rocblas_handle handle, const rocblas_fill uplo, const int m,
                                           Kokkos::complex<double> *a, const int lda, int *dev) {
  const int r_val = rocsolver_zpotrf(handle, uplo, m, (rocblas_double_complex *)a, lda, dev);
  return r_val;
}
#endif

template <>
int Lapack<Kokkos::complex<double>>::sytrf(const char uplo, const int m, Kokkos::complex<double> *a, const int lda,
                                           int *ipiv, Kokkos::complex<double> *work, int lwork, int *info) {
  F77_FUNC_ZSYTRF(&uplo, &m, a, &lda, ipiv, work, &lwork, info);
  return 0;
}

#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<Kokkos::complex<double>>::sytrf_buffersize(cusolverDnHandle_t handle, const int m,
                                                      Kokkos::complex<double> *a, const int lda, int *lwork) {
  const int r_val = cusolverDnZsytrf_bufferSize(handle, m, (cuDoubleComplex *)a, lda, lwork);
  return r_val;
}

template <>
int Lapack<Kokkos::complex<double>>::sytrf(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m,
                                           Kokkos::complex<double> *a, const int lda, int *ipiv,
                                           Kokkos::complex<double> *w, const int lwork, int *dev) {
  const int r_val =
      cusolverDnZsytrf(handle, uplo, m, (cuDoubleComplex *)a, lda, ipiv, (cuDoubleComplex *)w, lwork, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<Kokkos::complex<double>>::sytrf(rocblas_handle handle, const rocblas_fill uplo, const int m,
                                           Kokkos::complex<double> *a, const int lda, int *ipiv, int *dev) {
#if defined(TACHO_ENABLE_ROCSOLVER_SYTRF)
  const int r_val = rocsolver_zsytrf(handle, uplo, m, (rocblas_double_complex *)a, lda, ipiv, dev);
  return r_val;
#else
  throw std::logic_error("Error: sytrf is not avail in old ROCM");
  return -1;
#endif
}
#endif

template <>
int Lapack<Kokkos::complex<double>>::getrf(const int m, const int n, Kokkos::complex<double> *a, const int lda,
                                           int *ipiv, int *info) {
  F77_FUNC_ZGETRF(&m, &n, a, &lda, ipiv, info);
  return 0;
}

#if defined(TACHO_ENABLE_CUSOLVER)
template <>
int Lapack<Kokkos::complex<double>>::getrf_buffersize(cusolverDnHandle_t handle, const int m, const int n,
                                                      Kokkos::complex<double> *a, const int lda, int *lwork) {
  const int r_val = cusolverDnZgetrf_bufferSize(handle, m, n, (cuDoubleComplex *)a, lda, lwork);
  return r_val;
}

template <>
int Lapack<Kokkos::complex<double>>::getrf(cusolverDnHandle_t handle, const int m, const int n,
                                           Kokkos::complex<double> *a, const int lda, Kokkos::complex<double> *w,
                                           int *ipiv, int *dev) {
  const int r_val = cusolverDnZgetrf(handle, m, n, (cuDoubleComplex *)a, lda, (cuDoubleComplex *)w, ipiv, dev);
  return r_val;
}
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
template <>
int Lapack<Kokkos::complex<double>>::getrf(rocblas_handle handle, const int m, const int n, Kokkos::complex<double> *a,
                                           const int lda, int *ipiv, int *dev) {
  const int r_val = rocsolver_zgetrf(handle, m, n, (rocblas_double_complex *)a, lda, ipiv, dev);
  return r_val;
}
#endif

template <>
int Lapack<std::complex<float>>::potrf(const char uplo, const int m, std::complex<float> *a, const int lda, int *info) {
  F77_FUNC_CPOTRF(&uplo, &m, (Kokkos::complex<float> *)a, &lda, info);
  return 0;
}

template <>
int Lapack<std::complex<float>>::sytrf(const char uplo, const int m, std::complex<float> *a, const int lda, int *ipiv,
                                       std::complex<float> *work, int lwork, int *info) {
  F77_FUNC_CSYTRF(&uplo, &m, (Kokkos::complex<float> *)a, &lda, ipiv, (Kokkos::complex<float> *)work, &lwork, info);
  return 0;
}

template <>
int Lapack<std::complex<double>>::potrf(const char uplo, const int m, std::complex<double> *a, const int lda,
                                        int *info) {
  F77_FUNC_ZPOTRF(&uplo, &m, (Kokkos::complex<double> *)a, &lda, info);
  return 0;
}
template <>
int Lapack<std::complex<double>>::sytrf(const char uplo, const int m, std::complex<double> *a, const int lda, int *ipiv,
                                        std::complex<double> *work, int lwork, int *info) {
  F77_FUNC_ZSYTRF(&uplo, &m, (Kokkos::complex<double> *)a, &lda, ipiv, (Kokkos::complex<double> *)work, &lwork, info);
  return 0;
}

} // namespace Tacho
