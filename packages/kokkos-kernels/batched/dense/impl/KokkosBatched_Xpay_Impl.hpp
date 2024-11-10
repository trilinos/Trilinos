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
#ifndef __KOKKOSBATCHED_XPAY_IMPL_HPP__
#define __KOKKOSBATCHED_XPAY_IMPL_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
struct SerialXpayInternal {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const ScalarType alpha, const ValueType* KOKKOS_RESTRICT X,
                                           const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) {
      Y[i * ys0] *= alpha;
      Y[i * ys0] += X[i * xs0];
    }

    return 0;
  }

  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const ScalarType* KOKKOS_RESTRICT alpha, const int alphas0,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) {
      Y[i * ys0] *= alpha[i * alphas0];
      Y[i * ys0] += X[i * xs0];
    }

    return 0;
  }

  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int n, const ScalarType* KOKKOS_RESTRICT alpha,
                                           const int alphas0, const ValueType* KOKKOS_RESTRICT X, const int xs0,
                                           const int xs1,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0, const int ys1) {
    if (xs0 > xs1)
      for (int i = 0; i < m; ++i) invoke(n, alpha[i * alphas0], X + i * xs0, xs1, Y + i * ys0, ys1);
    else
      for (int j = 0; j < n; ++j) invoke(m, alpha, alphas0, X + j * xs1, xs0, Y + j * ys1, ys0);

    return 0;
  }
};

///
/// Team Internal Impl
/// ====================
struct TeamXpayInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const int m, const ScalarType alpha,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int& i) {
      Y[i * ys0] *= alpha;
      Y[i * ys0] += X[i * xs0];
    });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const int m,
                                           const ScalarType* KOKKOS_RESTRICT alpha, const int alphas0,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int& i) {
      Y[i * ys0] *= alpha[i * alphas0];
      Y[i * ys0] += X[i * xs0];
    });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const int m, const int n,
                                           const ScalarType* KOKKOS_RESTRICT alpha, const int alphas0,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0, const int xs1,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0, const int ys1) {
    if (m > n) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int& i) {
        SerialXpayInternal::invoke(n, alpha[i * alphas0], X + i * xs0, xs1, Y + i * ys0, ys1);
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int& j) {
        SerialXpayInternal::invoke(m, alpha, alphas0, X + j * xs1, xs0, Y + j * ys1, ys0);
      });
    }
    // member.team_barrier();
    return 0;
  }
};

///
/// TeamVector Internal Impl
/// ========================
struct TeamVectorXpayInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const int m, const ScalarType alpha,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int& i) {
      Y[i * ys0] *= alpha;
      Y[i * ys0] += X[i * xs0];
    });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const int m,
                                           const ScalarType* KOKKOS_RESTRICT alpha, const int alphas0,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int& i) {
      Y[i * ys0] *= alpha[i * alphas0];
      Y[i * ys0] += X[i * xs0];
    });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType, typename layout>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const int m, const int n,
                                           const ScalarType* KOKKOS_RESTRICT alpha, const int alphas0,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0, const int xs1,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0, const int ys1) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, m * n), [&](const int& iTemp) {
      int i, j;
      getIndices<int, layout>(iTemp, n, m, j, i);
      Y[i * ys0 + j * ys1] *= alpha[i * alphas0];
      Y[i * ys0 + j * ys1] += X[i * xs0 + j * xs1];
    });
    // member.team_barrier();
    return 0;
  }
};

///
/// Serial Impl
/// ===========
template <typename ViewType, typename alphaViewType>
KOKKOS_INLINE_FUNCTION int SerialXpay::invoke(const alphaViewType& alpha, const ViewType& X, const ViewType& Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert(Kokkos::is_view<ViewType>::value, "KokkosBatched::xpay: ViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<alphaViewType>::value, "KokkosBatched::xpay: alphaViewType is not a Kokkos::View.");
  static_assert(ViewType::rank == 2, "KokkosBatched::xpay: ViewType must have rank 2.");
  static_assert(alphaViewType::rank == 1, "KokkosBatched::xpay: alphaViewType must have rank 1.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
    Kokkos::printf(
        "KokkosBatched::xpay: Dimensions of X and Y do not match: X: %d x %d, "
        "Y: %d x %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
    return 1;
  }
  if (X.extent(0) != alpha.extent(0)) {
    Kokkos::printf(
        "KokkosBatched::xpay: First dimension of X and alpha do not match: X: "
        "%d x %d, alpha: %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)alpha.extent(0));
    return 1;
  }
#endif

  return SerialXpayInternal::template invoke<typename alphaViewType::non_const_value_type,
                                             typename ViewType::non_const_value_type>(
      X.extent(0), X.extent(1), alpha.data(), alpha.stride_0(), X.data(), X.stride_0(), X.stride_1(), Y.data(),
      Y.stride_0(), Y.stride_1());
}

///
/// Team Impl
/// =========

template <typename MemberType>
template <typename ViewType, typename alphaViewType>
KOKKOS_INLINE_FUNCTION int TeamXpay<MemberType>::invoke(const MemberType& member, const alphaViewType& alpha,
                                                        const ViewType& X, const ViewType& Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert(Kokkos::is_view<ViewType>::value, "KokkosBatched::xpay: ViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<alphaViewType>::value, "KokkosBatched::xpay: alphaViewType is not a Kokkos::View.");
  static_assert(ViewType::rank == 2, "KokkosBatched::xpay: ViewType must have rank 2.");
  static_assert(alphaViewType::rank == 1, "KokkosBatched::xpay: alphaViewType must have rank 1.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
    Kokkos::printf(
        "KokkosBatched::xpay: Dimensions of X and Y do not match: X: %d x %d, "
        "Y: %d x %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
    return 1;
  }
  if (X.extent(0) != alpha.extent(0)) {
    Kokkos::printf(
        "KokkosBatched::xpay: First dimension of X and alpha do not match: X: "
        "%d x %d, alpha: %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)alpha.extent(0));
    return 1;
  }
#endif

  return TeamXpayInternal::template invoke<MemberType, typename alphaViewType::non_const_value_type,
                                           typename ViewType::non_const_value_type>(
      member, X.extent(0), X.extent(1), alpha.data(), alpha.stride_0(), X.data(), X.stride_0(), X.stride_1(), Y.data(),
      Y.stride_0(), Y.stride_1());
}

///
/// TeamVector Impl
/// ===============

template <typename MemberType>
template <typename ViewType, typename alphaViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorXpay<MemberType>::invoke(const MemberType& member, const alphaViewType& alpha,
                                                              const ViewType& X, const ViewType& Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert(Kokkos::is_view<ViewType>::value, "KokkosBatched::xpay: ViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<alphaViewType>::value, "KokkosBatched::xpay: alphaViewType is not a Kokkos::View.");
  static_assert(ViewType::rank == 2, "KokkosBatched::xpay: ViewType must have rank 2.");
  static_assert(alphaViewType::rank == 1, "KokkosBatched::xpay: alphaViewType must have rank 1.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
    Kokkos::printf(
        "KokkosBatched::xpay: Dimensions of X and Y do not match: X: %d x %d, "
        "Y: %d x %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
    return 1;
  }
  if (X.extent(0) != alpha.extent(0)) {
    Kokkos::printf(
        "KokkosBatched::xpay: First dimension of X and alpha do not match: X: "
        "%d x %d, alpha: %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)alpha.extent(0));
    return 1;
  }
#endif

  return TeamVectorXpayInternal::invoke<MemberType, typename alphaViewType::non_const_value_type,
                                        typename ViewType::non_const_value_type, typename ViewType::array_layout>(
      member, X.extent(0), X.extent(1), alpha.data(), alpha.stride_0(), X.data(), X.stride_0(), X.stride_1(), Y.data(),
      Y.stride_0(), Y.stride_1());
}

}  // namespace KokkosBatched

#endif
