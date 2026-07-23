// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_ROTM_INTERNAL_HPP_
#define KOKKOSBATCHED_ROTM_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================
template <int Flag>
struct SerialRotmInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static void invoke(const int n, ValueType *KOKKOS_RESTRICT x, const int xs0,
                                            ValueType *KOKKOS_RESTRICT y, const int ys0,
                                            const ValueType *KOKKOS_RESTRICT param, const int ps0);
};

template <int Flag>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION void SerialRotmInternal<Flag>::invoke(const int n, ValueType *KOKKOS_RESTRICT x, const int xs0,
                                                             ValueType *KOKKOS_RESTRICT y, const int ys0,
                                                             const ValueType *KOKKOS_RESTRICT param, const int ps0) {
  if constexpr (Flag == -1) {
    // flag == -1.0: [[h11, h12], [h21, h22]]
    auto h11 = param[1 * ps0];
    auto h21 = param[2 * ps0];
    auto h12 = param[3 * ps0];
    auto h22 = param[4 * ps0];
    for (int i = 0; i < n; ++i) {
      auto temp  = h11 * x[i * xs0] + h12 * y[i * ys0];
      y[i * ys0] = h21 * x[i * xs0] + h22 * y[i * ys0];
      x[i * xs0] = temp;
    }
  } else if constexpr (Flag == 0) {
    // flag == 0.0: [[1, h12], [h21, 1]]
    auto h12 = param[3 * ps0];
    auto h21 = param[2 * ps0];
    for (int i = 0; i < n; ++i) {
      auto temp  = x[i * xs0] + h12 * y[i * ys0];
      y[i * ys0] = h21 * x[i * xs0] + y[i * ys0];
      x[i * xs0] = temp;
    }
  } else if constexpr (Flag == 1) {
    // flag == 1.0: [[h11, 1], [-1, h22]]
    auto h11 = param[1 * ps0];
    auto h22 = param[4 * ps0];
    for (int i = 0; i < n; ++i) {
      auto temp  = h11 * x[i * xs0] + y[i * ys0];
      y[i * ys0] = -x[i * xs0] + h22 * y[i * ys0];
      x[i * xs0] = temp;
    }
  }
  // Checks for flag == -2.0 (identity) are done in the checkRotmInput function, so we don't need to handle it here.
}

///
/// Team Internal Impl
///
/// ==================
template <int Flag>
struct TeamRotmInternal {
  template <typename MemberType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static void invoke(const MemberType &member, const int n, ValueType *KOKKOS_RESTRICT x,
                                            const int xs0, ValueType *KOKKOS_RESTRICT y, const int ys0,
                                            const ValueType *KOKKOS_RESTRICT param, const int ps0);
};

template <int Flag>
template <typename MemberType, typename ValueType>
KOKKOS_INLINE_FUNCTION void TeamRotmInternal<Flag>::invoke(const MemberType &member, const int n,
                                                           ValueType *KOKKOS_RESTRICT x, const int xs0,
                                                           ValueType *KOKKOS_RESTRICT y, const int ys0,
                                                           const ValueType *KOKKOS_RESTRICT param, const int ps0) {
  if constexpr (Flag == -1) {
    // flag == -1.0: [[h11, h12], [h21, h22]]
    auto h11 = param[1 * ps0];
    auto h21 = param[2 * ps0];
    auto h12 = param[3 * ps0];
    auto h22 = param[4 * ps0];
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &i) {
      auto temp  = h11 * x[i * xs0] + h12 * y[i * ys0];
      y[i * ys0] = h21 * x[i * xs0] + h22 * y[i * ys0];
      x[i * xs0] = temp;
    });
  } else if constexpr (Flag == 0) {
    // flag == 0.0: [[1, h12], [h21, 1]]
    auto h12 = param[3 * ps0];
    auto h21 = param[2 * ps0];
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &i) {
      auto temp  = x[i * xs0] + h12 * y[i * ys0];
      y[i * ys0] = h21 * x[i * xs0] + y[i * ys0];
      x[i * xs0] = temp;
    });
  } else if constexpr (Flag == 1) {
    // flag == 1.0: [[h11, 1], [-1, h22]]
    auto h11 = param[1 * ps0];
    auto h22 = param[4 * ps0];
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &i) {
      auto temp  = h11 * x[i * xs0] + y[i * ys0];
      y[i * ys0] = -x[i * xs0] + h22 * y[i * ys0];
      x[i * xs0] = temp;
    });
  }
}

///
/// TeamVector Internal Impl
/// ========================
template <int Flag>
struct TeamVectorRotmInternal {
  template <typename MemberType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static void invoke(const MemberType &member, const int n, ValueType *KOKKOS_RESTRICT x,
                                            const int xs0, ValueType *KOKKOS_RESTRICT y, const int ys0,
                                            const ValueType *KOKKOS_RESTRICT param, const int ps0);
};

template <int Flag>
template <typename MemberType, typename ValueType>
KOKKOS_INLINE_FUNCTION void TeamVectorRotmInternal<Flag>::invoke(const MemberType &member, const int n,
                                                                 ValueType *KOKKOS_RESTRICT x, const int xs0,
                                                                 ValueType *KOKKOS_RESTRICT y, const int ys0,
                                                                 const ValueType *KOKKOS_RESTRICT param,
                                                                 const int ps0) {
  if constexpr (Flag == -1) {
    // flag == -1.0: [[h11, h12], [h21, h22]]
    auto h11 = param[1 * ps0];
    auto h21 = param[2 * ps0];
    auto h12 = param[3 * ps0];
    auto h22 = param[4 * ps0];
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int &i) {
      auto temp  = h11 * x[i * xs0] + h12 * y[i * ys0];
      y[i * ys0] = h21 * x[i * xs0] + h22 * y[i * ys0];
      x[i * xs0] = temp;
    });
  } else if constexpr (Flag == 0) {
    // flag == 0.0: [[1, h12], [h21, 1]]
    auto h12 = param[3 * ps0];
    auto h21 = param[2 * ps0];
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int &i) {
      auto temp  = x[i * xs0] + h12 * y[i * ys0];
      y[i * ys0] = h21 * x[i * xs0] + y[i * ys0];
      x[i * xs0] = temp;
    });
  } else if constexpr (Flag == 1) {
    // flag == 1.0: [[h11, 1], [-1, h22]]
    auto h11 = param[1 * ps0];
    auto h22 = param[4 * ps0];
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n), [&](const int &i) {
      auto temp  = h11 * x[i * xs0] + y[i * ys0];
      y[i * ys0] = -x[i * xs0] + h22 * y[i * ys0];
      x[i * xs0] = temp;
    });
  }
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_ROTM_INTERNAL_HPP_
