// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_TRMV_INTERNAL_HPP_
#define KOKKOSBATCHED_TRMV_INTERNAL_HPP_

#include <concepts>
#include <KokkosBatched_Util.hpp>

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

struct SerialTrmvInternal {
  template <typename ArgUplo, typename ArgDiag, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int n, const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
                                           const int as0, const int as1, ValueType *KOKKOS_RESTRICT x, const int xs0);

  template <typename ArgUplo, typename ArgDiag, typename Op, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, const int n, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           ValueType *KOKKOS_RESTRICT x, const int xs0);
};

// Non-transpose case
template <typename ArgUplo, typename ArgDiag, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrmvInternal::invoke(const int n, const ScalarType alpha,
                                                      const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                                      ValueType *KOKKOS_RESTRICT x, const int xs0) {
  if constexpr (std::same_as<ArgUplo, KokkosBatched::Uplo::Lower>) {
    // Lower triangular specific implementation
    for (int j = n - 1; j >= 0; j--) {
      if (x[j * xs0] != ScalarType(0)) {
        x[j * xs0] *= alpha;
        auto temp = x[j * xs0];
        for (int i = n - 1; i > j; i--) {
          x[i * xs0] += temp * A[i * as0 + j * as1];
        }
        if constexpr (std::same_as<ArgDiag, KokkosBatched::Diag::NonUnit>) {
          x[j * xs0] *= A[j * as0 + j * as1];
        }
      }
    }
  } else {
    // Upper triangular specific implementation
    for (int j = 0; j < n; j++) {
      if (x[j * xs0] != ScalarType(0)) {
        x[j * xs0] *= alpha;
        auto temp = x[j * xs0];
        for (int i = 0; i < j; i++) {
          x[i * xs0] += temp * A[i * as0 + j * as1];
        }
        if constexpr (std::same_as<ArgDiag, KokkosBatched::Diag::NonUnit>) {
          x[j * xs0] *= A[j * as0 + j * as1];
        }
      }
    }
  }

  return 0;
}

// Transpose or Conjugate-transpose case
template <typename ArgUplo, typename ArgDiag, typename Op, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTrmvInternal::invoke(Op op, const int n, const ScalarType alpha,
                                                      const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                                      ValueType *KOKKOS_RESTRICT x, const int xs0) {
  if constexpr (std::same_as<ArgUplo, KokkosBatched::Uplo::Lower>) {
    // Lower triangular specific implementation
    for (int j = 0; j < n; j++) {
      auto temp = x[j * xs0];
      if constexpr (std::same_as<ArgDiag, KokkosBatched::Diag::NonUnit>) {
        temp *= op(A[j * as0 + j * as1]);
      }
      for (int i = j + 1; i < n; i++) {
        temp += op(A[i * as0 + j * as1]) * x[i * xs0];
      }
      x[j * xs0] = temp * alpha;
    }
  } else {
    // Upper triangular specific implementation
    for (int j = n - 1; j >= 0; j--) {
      auto temp = x[j * xs0];
      if constexpr (std::same_as<ArgDiag, KokkosBatched::Diag::NonUnit>) {
        temp *= op(A[j * as0 + j * as1]);
      }
      for (int i = j - 1; i >= 0; i--) {
        temp += op(A[i * as0 + j * as1]) * x[i * xs0];
      }
      x[j * xs0] = temp * alpha;
    }
  }

  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_TRMV_INTERNAL_HPP_
