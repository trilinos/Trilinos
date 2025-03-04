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

#ifndef KOKKOSBATCHED_GER_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_GER_SERIAL_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

struct SerialGerInternal {
  template <typename Op, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, const int am, const int an, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT x, const int xs0,
                                           const ValueType *KOKKOS_RESTRICT y, const int ys0,
                                           ValueType *KOKKOS_RESTRICT A, const int as0, const int as1);
};

template <typename Op, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialGerInternal::invoke(Op op, const int am, const int an, const ScalarType alpha,
                                                     const ValueType *KOKKOS_RESTRICT x, const int xs0,
                                                     const ValueType *KOKKOS_RESTRICT y, const int ys0,
                                                     ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
  for (int j = 0; j < an; j++) {
    if (y[j * ys0] != 0) {
      auto temp = alpha * op(y[j * ys0]);
      for (int i = 0; i < am; i++) {
        A[i * as0 + j * as1] += x[i * xs0] * temp;
      }
    }
  }

  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_GER_SERIAL_INTERNAL_HPP_
