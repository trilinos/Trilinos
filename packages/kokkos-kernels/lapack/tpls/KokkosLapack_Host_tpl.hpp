// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_HOST_TPL_HPP_
#define KOKKOSLAPACK_HOST_TPL_HPP_

/// \file  KokkosLapack_Host_tpl.hpp
/// \brief LAPACK wrapper

#include "KokkosKernels_config.h"
#include "KokkosKernels_ArithTraits.hpp"

#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK)

namespace KokkosLapack {
namespace Impl {

template <typename T>
struct HostLapack {
  static void gesv(int n, int rhs, T *a, int lda, int *ipiv, T *b, int ldb, int info);

  static void gesvd(const char jobu, const char jobvt, const int m, const int n, T *A, const int lda,
                    typename KokkosKernels::ArithTraits<T>::mag_type *S, T *U, const int ldu, T *Vt, const int ldvt,
                    T *work, int lwork, typename KokkosKernels::ArithTraits<T>::mag_type *rwork, int info);

  static int trtri(const char uplo, const char diag, int n, const T *a, int lda);
};
}  // namespace Impl
}  // namespace KokkosLapack

#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK

#endif  // KOKKOSLAPACK_HOST_TPL_HPP_
