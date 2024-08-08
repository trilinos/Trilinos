//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSBLAS_HOST_TPL_HPP_
#define KOKKOSBLAS_HOST_TPL_HPP_

/// \file  KokkosBlas_Host_tpl.hpp
/// \brief BLAS wrapper
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_config.h"
#include "Kokkos_ArithTraits.hpp"

#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)
#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
#include "mkl_types.h"
#endif

namespace KokkosBlas {
namespace Impl {

#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
using KK_INT = MKL_INT;
#else
using KK_INT = int;
#endif

template <typename T>
struct HostBlas {
  typedef Kokkos::ArithTraits<T> ats;
  typedef typename ats::mag_type mag_type;

  static void scal(KK_INT n, const T alpha,
                   /* */ T *x, KK_INT x_inc);

  static KK_INT iamax(KK_INT n, const T *x, KK_INT x_inc);

  static mag_type nrm2(KK_INT n, const T *x, KK_INT x_inc);

  static mag_type asum(KK_INT n, const T *x, KK_INT x_inc);

  static T dot(KK_INT n, const T *x, KK_INT x_inc, const T *y, KK_INT y_inc);

  static void axpy(KK_INT n, const T alpha, const T *x, KK_INT x_inc,
                   /* */ T *y, KK_INT y_inc);

  static void rot(KK_INT const N, T *X, KK_INT const incx, T *Y, KK_INT const incy, mag_type *c, mag_type *s);

  static void rotg(T *a, T *b, mag_type *c, T *s);

  static void rotm(const KK_INT n, T *X, const KK_INT incx, T *Y, const KK_INT incy, T const *param);

  static void rotmg(T *d1, T *d2, T *x1, const T *y1, T *param);

  static void swap(KK_INT const N, T *X, KK_INT const incx, T *Y, KK_INT const incy);

  static void gemv(const char trans, KK_INT m, KK_INT n, const T alpha, const T *a, KK_INT lda, const T *b, KK_INT ldb,
                   const T beta,
                   /* */ T *c, KK_INT ldc);

  static void ger(KK_INT m, KK_INT n, const T alpha, const T *x, KK_INT incx, const T *y, KK_INT incy, T *a,
                  KK_INT lda);

  static void geru(KK_INT m, KK_INT n, const T alpha, const T *x, KK_INT incx, const T *y, KK_INT incy, T *a,
                   KK_INT lda);

  static void gerc(KK_INT m, KK_INT n, const T alpha, const T *x, KK_INT incx, const T *y, KK_INT incy, T *a,
                   KK_INT lda);

  static void syr(const char uplo, KK_INT n, const T alpha, const T *x, KK_INT incx, T *a, KK_INT lda);

  static void syr2(const char uplo, KK_INT n, const T alpha, const T *x, KK_INT incx, const T *y, KK_INT incy, T *a,
                   KK_INT lda);

  template <typename tAlpha>
  static void her(const char uplo, KK_INT n, const tAlpha alpha, const T *x, KK_INT incx, T *a, KK_INT lda);

  static void her2(const char uplo, KK_INT n, const T alpha, const T *x, KK_INT incx, const T *y, KK_INT incy, T *a,
                   KK_INT lda);

  static void trsv(const char uplo, const char transa, const char diag, KK_INT m, const T *a, KK_INT lda,
                   /* */ T *b, KK_INT ldb);

  static void gemm(const char transa, const char transb, KK_INT m, KK_INT n, KK_INT k, const T alpha, const T *a,
                   KK_INT lda, const T *b, KK_INT ldb, const T beta,
                   /* */ T *c, KK_INT ldc);

  static void herk(const char transa, const char transb, KK_INT n, KK_INT k, const T alpha, const T *a, KK_INT lda,
                   const T beta,
                   /* */ T *c, KK_INT ldc);

  static void trmm(const char side, const char uplo, const char transa, const char diag, KK_INT m, KK_INT n,
                   const T alpha, const T *a, KK_INT lda,
                   /* */ T *b, KK_INT ldb);

  static void trsm(const char side, const char uplo, const char transa, const char diag, KK_INT m, KK_INT n,
                   const T alpha, const T *a, KK_INT lda,
                   /* */ T *b, KK_INT ldb);
};
}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

#endif  // KOKKOSBLAS_HOST_TPL_HPP_
