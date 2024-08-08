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

#ifndef KOKKOSLAPACK_HOST_TPL_HPP_
#define KOKKOSLAPACK_HOST_TPL_HPP_

/// \file  KokkosLapack_Host_tpl.hpp
/// \brief LAPACK wrapper

#include "KokkosKernels_config.h"
#include "Kokkos_ArithTraits.hpp"

#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK)

namespace KokkosLapack {
namespace Impl {

template <typename T>
struct HostLapack {
  static void gesv(int n, int rhs, T *a, int lda, int *ipiv, T *b, int ldb, int info);

  static void gesvd(const char jobu, const char jobvt, const int m, const int n, T *A, const int lda,
                    typename Kokkos::ArithTraits<T>::mag_type *S, T *U, const int ldu, T *Vt, const int ldvt, T *work,
                    int lwork, typename Kokkos::ArithTraits<T>::mag_type *rwork, int info);

  static int trtri(const char uplo, const char diag, int n, const T *a, int lda);
};
}  // namespace Impl
}  // namespace KokkosLapack

#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK

#endif  // KOKKOSLAPACK_HOST_TPL_HPP_
