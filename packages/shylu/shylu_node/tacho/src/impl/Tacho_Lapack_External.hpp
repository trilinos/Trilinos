#ifndef __TACHO_LAPACK_EXTERNAL_HPP__
#define __TACHO_LAPACK_EXTERNAL_HPP__

/// \file  Tacho_Blas_External.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Core.hpp"

#if defined(KOKKOS_ENABLE_CUDA)
#define TACHO_ENABLE_CUSOLVER
#endif

#if defined(KOKKOS_ENABLE_HIP)
// todo: enable hipblas interface after checking on AMD machine
//#define TACHO_ENABLE_HIPSOLVER
#endif

#if defined(TACHO_ENABLE_CUSOLVER)
#include "cusolverDn.h"
#endif
#if defined(TACHO_ENABLE_HIPSOLVER)
#include "hipsolver.h"
#endif

namespace Tacho {

template <typename T> struct Lapack {
  ///
  /// Cholesky
  ///
  static int potrf(const char uplo, const int m, T *a, const int lda, int *info);
#if defined(TACHO_ENABLE_CUSOLVER)
  static int potrf_buffersize(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m, T *a, const int lda,
                              int *lwork);
  static int potrf(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m, T *a, const int lda, T *W,
                   const int lwork, int *dev);
#endif
#if defined(TACHO_ENABLE_HIPSOLVER)
  static int potrf_buffersize(hipsolverHandle_t handle, const hipblasFillMode_t uplo, const int m, T *a, const int lda,
                              int *lwork);
  static int potrf(hipsolverHandle_t handle, const hipblasFillMode_t uplo, const int m, T *a, const int lda, T *W,
                   const int lwork, int *dev);
#endif

  ///
  /// LDLt
  ///
  static int sytrf(const char uplo, const int m, T *a, const int lda, int *ipiv, T *work, int lwork, int *info);
#if defined(TACHO_ENABLE_CUSOLVER)
  static int sytrf_buffersize(cusolverDnHandle_t handle, const int m, T *a, const int lda, int *lwork);
  static int sytrf(cusolverDnHandle_t handle, const cublasFillMode_t uplo, const int m, T *a, const int lda, int *ipiv,
                   T *W, const int lwork, int *dev);
#endif
#if defined(TACHO_ENABLE_HIPSOLVER)
  static int sytrf_buffersize(hipsolverHandle_t handle, const int m, T *a, const int lda, int *lwork);
  static int sytrf(hipsolverHandle_t handle, const hipblasFillMode_t uplo, const int m, T *a, const int lda, int *ipiv,
                   T *W, const int lwork, int *dev);
#endif

  ///
  /// LU
  ///
  static int getrf(const int m, const int n, T *a, const int lda, int *ipiv, int *info);
#if defined(TACHO_ENABLE_CUSOLVER)
  static int getrf_buffersize(cusolverDnHandle_t handle, const int m, const int n, T *a, const int lda, int *lwork);
  static int getrf(cusolverDnHandle_t handle, const int m, const int n, T *a, const int lda, T *w, int *ipiv, int *dev);
#endif
#if defined(TACHO_ENABLE_HIPSOLVER)
  static int getrf_buffersize(hipsolverHandle_t handle, const int m, const int n, T *a, const int lda, int *lwork);
  static int getrf(hipsolverHandle_t handle, const int m, const int n, T *a, const int lda, T *w, int *ipiv, int *dev);
#endif
};
} // namespace Tacho

#endif
