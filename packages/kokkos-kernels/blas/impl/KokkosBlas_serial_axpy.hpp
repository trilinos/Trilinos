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

#ifndef KOKKOSBLAS_SERIAL_AXPY_IMPL_HPP_
#define KOKKOSBLAS_SERIAL_AXPY_IMPL_HPP_

#include <Kokkos_Core.hpp>

namespace KokkosBlas {
namespace Impl {

///
/// Serial Internal Impl
/// ====================
template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION static void serial_axpy(const int m, const ScalarType alpha, const ValueType *KOKKOS_RESTRICT X,
                                               /* */ ValueType *KOKKOS_RESTRICT Y, const int xs0, const int ys0) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int i = 0; i < m; ++i) Y[i * ys0] += alpha * X[i * xs0];

  return;
}

template <typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION static void serial_axpy_mv(const int m, const int n, const ScalarType alpha,
                                                  const ValueType *KOKKOS_RESTRICT X,
                                                  /* */ ValueType *KOKKOS_RESTRICT Y, const int xs0, const int xs1,
                                                  const int ys0, const int ys1) {
  if (xs0 > xs1) {
    for (int i = 0; i < m; ++i) serial_axpy(n, alpha, X + i * xs0, Y + i * ys0, xs1, ys1);
  } else {
    for (int j = 0; j < n; ++j) serial_axpy(m, alpha, X + j * xs1, Y + j * ys1, xs0, ys0);
  }

  return;
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif
