// clang-format off
/* =====================================================================================
Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

SCR#:2790.0

This file is part of Tacho. Tacho is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Kyungjoo Kim at <kyukim@sandia.gov,https://github.com/kyungjoo-kim>

Sandia National Laboratories, Albuquerque, NM, USA
===================================================================================== */
// clang-format on
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
#define TACHO_ENABLE_ROCSOLVER
#endif

#if defined(TACHO_ENABLE_CUSOLVER)
#include "cusolverDn.h"
#endif

#if defined(TACHO_ENABLE_ROCSOLVER)
#include "rocblas/rocblas.h"
#include "rocsolver/rocsolver.h"
#if (ROCSOLVER_VERSION_MAJOR >= 3) && (ROCSOLVER_VERSION_MINOR >= 16)
#define TACHO_ENABLE_ROCSOLVER_SYTRF 
#endif
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
#if defined(TACHO_ENABLE_ROCSOLVER)
  static int potrf(rocblas_handle handle, const rocblas_fill uplo, const int m, T *a, const int lda, int *dev);
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
#if defined(TACHO_ENABLE_ROCSOLVER)
  static int sytrf(rocblas_handle handle, const rocblas_fill uplo, const int m, T *a, const int lda, int *ipiv,
                   int *dev);
#endif

  ///
  /// LU
  ///
  static int getrf(const int m, const int n, T *a, const int lda, int *ipiv, int *info);
#if defined(TACHO_ENABLE_CUSOLVER)
  static int getrf_buffersize(cusolverDnHandle_t handle, const int m, const int n, T *a, const int lda, int *lwork);
  static int getrf(cusolverDnHandle_t handle, const int m, const int n, T *a, const int lda, T *w, int *ipiv, int *dev);
#endif
#if defined(TACHO_ENABLE_ROCSOLVER)
  static int getrf(rocblas_handle handle, const int m, const int n, T *a, const int lda, int *ipiv, int *dev);
#endif
};
} // namespace Tacho

#endif
