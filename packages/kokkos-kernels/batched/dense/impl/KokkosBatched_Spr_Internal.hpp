// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_SPR_INTERNAL_HPP_
#define KOKKOSBATCHED_SPR_INTERNAL_HPP_

#include <concepts>
#include <KokkosBatched_Util.hpp>

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

struct SerialSprInternal {
  template <typename ArgUplo, typename Op, typename SymOp, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, SymOp sym_op, const int n, const ScalarType alpha,
                                           const ValueType *KOKKOS_RESTRICT x, const int xs0,
                                           ValueType *KOKKOS_RESTRICT A, const int as0);
};

template <typename ArgUplo, typename Op, typename SymOp, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialSprInternal::invoke(Op op, SymOp sym_op, const int n, const ScalarType alpha,
                                                     const ValueType *KOKKOS_RESTRICT x, const int xs0,
                                                     ValueType *KOKKOS_RESTRICT ap, const int as0) {
  if constexpr (std::same_as<ArgUplo, Uplo::Lower>) {
    // Lower
    for (int j = 0; j < n; j++) {
      const int offset = j * (2 * n - j + 1) / 2;
      if (x[j * xs0] != ValueType(0)) {
        auto temp        = alpha * op(x[j * xs0]);
        ap[offset * as0] = sym_op(ap[offset * as0] + x[j * xs0] * temp);
        for (int i = j + 1; i < n; i++) {
          ap[(offset + i - j) * as0] += x[i * xs0] * temp;
        }
      } else {
        ap[offset * as0] = sym_op(ap[offset * as0]);
      }
    }
  } else {
    // Upper
    for (int j = 0; j < n; j++) {
      const int offset = j * (j + 1) / 2;
      if (x[j * xs0] != ValueType(0)) {
        auto temp = alpha * op(x[j * xs0]);
        for (int i = 0; i < j; i++) {
          ap[(offset + i) * as0] += x[i * xs0] * temp;
        }
        ap[(offset + j) * as0] = sym_op(ap[(offset + j) * as0] + x[j * xs0] * temp);
      } else {
        ap[(offset + j) * as0] = sym_op(ap[(offset + j) * as0]);
      }
    }
  }

  return 0;
}

///
/// Team Internal Impl
/// ====================

struct TeamSprInternal {
  template <typename ArgUplo, typename MemberType, typename Op, typename SymOp, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int n,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT x, const int xs0,
                                           ValueType *KOKKOS_RESTRICT A, const int as0);
};

template <typename ArgUplo, typename MemberType, typename Op, typename SymOp, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamSprInternal::invoke(const MemberType &member, Op op, SymOp sym_op, const int n,
                                                   const ScalarType alpha, const ValueType *KOKKOS_RESTRICT x,
                                                   const int xs0, ValueType *KOKKOS_RESTRICT ap, const int as0) {
  if constexpr (std::same_as<ArgUplo, Uplo::Lower>) {
    // Lower
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
      const int offset = j * (2 * n - j + 1) / 2;
      if (x[j * xs0] != ValueType(0)) {
        auto temp        = alpha * op(x[j * xs0]);
        ap[offset * as0] = sym_op(ap[offset * as0] + x[j * xs0] * temp);
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, (n - j - 1)),
                             [&](const int i) { ap[(offset + i + 1) * as0] += x[(i + j + 1) * xs0] * temp; });
      } else {
        ap[offset * as0] = sym_op(ap[offset * as0]);
      }
    });
  } else {
    // Upper
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int j) {
      const int offset = j * (j + 1) / 2;
      if (x[j * xs0] != ValueType(0)) {
        auto temp = alpha * op(x[j * xs0]);
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, j),
                             [&](const int i) { ap[(offset + i) * as0] += x[i * xs0] * temp; });
        ap[(offset + j) * as0] = sym_op(ap[(offset + j) * as0] + x[j * xs0] * temp);
      } else {
        ap[(offset + j) * as0] = sym_op(ap[(offset + j) * as0]);
      }
    });
  }

  return 0;
}

///
/// TeamVector Internal Impl
/// ====================

struct TeamVectorSprInternal {
  template <typename ArgUplo, typename MemberType, typename Op, typename SymOp, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, SymOp sym_op, const int n,
                                           const ScalarType alpha, const ValueType *KOKKOS_RESTRICT x, const int xs0,
                                           ValueType *KOKKOS_RESTRICT A, const int as0);
};

template <typename ArgUplo, typename MemberType, typename Op, typename SymOp, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamVectorSprInternal::invoke(const MemberType &member, Op op, SymOp sym_op, const int n,
                                                         const ScalarType alpha, const ValueType *KOKKOS_RESTRICT x,
                                                         const int xs0, ValueType *KOKKOS_RESTRICT ap, const int as0) {
  if constexpr (std::same_as<ArgUplo, Uplo::Lower>) {
    // Lower
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
      const int offset = j * (2 * n - j + 1) / 2;
      if (x[j * xs0] != ValueType(0)) {
        auto temp        = alpha * op(x[j * xs0]);
        ap[offset * as0] = sym_op(ap[offset * as0] + x[j * xs0] * temp);
        for (int i = j + 1; i < n; i++) {
          ap[(offset + i - j) * as0] += x[i * xs0] * temp;
        }
      } else {
        ap[offset * as0] = sym_op(ap[offset * as0]);
      }
    });
  } else {
    // Upper
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int j) {
      const int offset = j * (j + 1) / 2;
      if (x[j * xs0] != ValueType(0)) {
        auto temp = alpha * op(x[j * xs0]);
        for (int i = 0; i < j; i++) {
          ap[(offset + i) * as0] += x[i * xs0] * temp;
        }
        ap[(offset + j) * as0] = sym_op(ap[(offset + j) * as0] + x[j * xs0] * temp);
      } else {
        ap[(offset + j) * as0] = sym_op(ap[(offset + j) * as0]);
      }
    });
  }

  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_SPR_INTERNAL_HPP_
