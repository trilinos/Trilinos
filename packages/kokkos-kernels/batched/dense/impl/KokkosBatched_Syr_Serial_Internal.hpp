// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_SYR_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_SYR_SERIAL_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

/// Lower

struct SerialSyrInternalLower {
  template <typename Op, typename SymOp, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, SymOp sym_op, const int an, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT x, const int xs0,
                                           ValueType *KOKKOS_RESTRICT A, const int as0, const int as1);
};

template <typename Op, typename SymOp, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialSyrInternalLower::invoke(Op op, SymOp sym_op, const int an, const ScalarType alpha,
                                                          const ValueType *KOKKOS_RESTRICT x, const int xs0,
                                                          ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
  for (int j = 0; j < an; j++) {
    if (x[j * xs0] != ValueType(0)) {
      auto temp            = alpha * op(x[j * xs0]);
      A[j * as0 + j * as1] = sym_op(A[j * as0 + j * as1] + x[j * xs0] * temp);
      for (int i = j + 1; i < an; i++) {
        A[i * as0 + j * as1] += x[i * xs0] * temp;
      }
    } else {
      A[j * as0 + j * as1] = sym_op(A[j * as0 + j * as1]);
    }
  }

  return 0;
}

/// Upper

struct SerialSyrInternalUpper {
  template <typename Op, typename SymOp, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, SymOp sym_op, const int an, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT x, const int xs0,
                                           ValueType *KOKKOS_RESTRICT A, const int as0, const int as1);
};

template <typename Op, typename SymOp, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialSyrInternalUpper::invoke(Op op, SymOp sym_op, const int an, const ScalarType alpha,
                                                          const ValueType *KOKKOS_RESTRICT x, const int xs0,
                                                          ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
  for (int j = 0; j < an; j++) {
    if (x[j * xs0] != ValueType(0)) {
      auto temp = alpha * op(x[j * xs0]);
      for (int i = 0; i < j; i++) {
        A[i * as0 + j * as1] += x[i * xs0] * temp;
      }
      A[j * as0 + j * as1] = sym_op(A[j * as0 + j * as1] + x[j * xs0] * temp);
    } else {
      A[j * as0 + j * as1] = sym_op(A[j * as0 + j * as1]);
    }
  }

  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_SYR_SERIAL_INTERNAL_HPP_
