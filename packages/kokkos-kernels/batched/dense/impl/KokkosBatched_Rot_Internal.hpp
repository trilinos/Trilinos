// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_ROT_INTERNAL_HPP_
#define KOKKOSBATCHED_ROT_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

struct SerialRotInternal {
  template <typename Op, typename ValueType, typename CType, typename SType>
  KOKKOS_INLINE_FUNCTION static int invoke(Op op, const int n, ValueType *KOKKOS_RESTRICT x, const int xs0,
                                           ValueType *KOKKOS_RESTRICT y, const int ys0, const CType c, const SType s);
};

template <typename Op, typename ValueType, typename CType, typename SType>
KOKKOS_INLINE_FUNCTION int SerialRotInternal::invoke(Op op, const int n, ValueType *KOKKOS_RESTRICT x, const int xs0,
                                                     ValueType *KOKKOS_RESTRICT y, const int ys0, const CType c,
                                                     const SType s) {
  for (int i = 0; i < n; i++) {
    auto temp  = c * x[i * xs0] + s * y[i * ys0];
    y[i * ys0] = c * y[i * ys0] - op(s) * x[i * xs0];
    x[i * xs0] = temp;
  }

  return 0;
}

///
/// Team Internal Impl
/// ==================

struct TeamRotInternal {
  template <typename MemberType, typename Op, typename ValueType, typename CType, typename SType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, const int n, ValueType *KOKKOS_RESTRICT x,
                                           const int xs0, ValueType *KOKKOS_RESTRICT y, const int ys0, const CType c,
                                           const SType s);
};

template <typename MemberType, typename Op, typename ValueType, typename CType, typename SType>
KOKKOS_INLINE_FUNCTION int TeamRotInternal::invoke(const MemberType &member, Op op, const int n,
                                                   ValueType *KOKKOS_RESTRICT x, const int xs0,
                                                   ValueType *KOKKOS_RESTRICT y, const int ys0, const CType c,
                                                   const SType s) {
  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &i) {
    auto temp  = c * x[i * xs0] + s * y[i * ys0];
    y[i * ys0] = c * y[i * ys0] - op(s) * x[i * xs0];
    x[i * xs0] = temp;
  });
  return 0;
}

///
/// TeamVector Internal Impl
/// ========================

struct TeamVectorRotInternal {
  template <typename MemberType, typename Op, typename ValueType, typename CType, typename SType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, Op op, const int n, ValueType *KOKKOS_RESTRICT x,
                                           const int xs0, ValueType *KOKKOS_RESTRICT y, const int ys0, const CType c,
                                           const SType s);
};

template <typename MemberType, typename Op, typename ValueType, typename CType, typename SType>
KOKKOS_INLINE_FUNCTION int TeamVectorRotInternal::invoke(const MemberType &member, Op op, const int n,
                                                         ValueType *KOKKOS_RESTRICT x, const int xs0,
                                                         ValueType *KOKKOS_RESTRICT y, const int ys0, const CType c,
                                                         const SType s) {
  Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int &i) {
    auto temp  = c * x[i * xs0] + s * y[i * ys0];
    y[i * ys0] = c * y[i * ys0] - op(s) * x[i * xs0];
    x[i * xs0] = temp;
  });
  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_ROT_INTERNAL_HPP_
