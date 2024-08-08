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
#ifndef __KOKKOSBATCHED_AXPY_IMPL_HPP__
#define __KOKKOSBATCHED_AXPY_IMPL_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBlas1_team_axpby.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
struct SerialAxpyInternal {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const ScalarType alpha, const ValueType* KOKKOS_RESTRICT X,
                                           const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) Y[i * ys0] += alpha * X[i * xs0];

    return 0;
  }

  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const ScalarType* KOKKOS_RESTRICT alpha, const int alphas0,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) Y[i * ys0] += alpha[i * alphas0] * X[i * xs0];

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
struct TeamAxpyInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const int m, const ScalarType alpha,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int& i) { Y[i * ys0] += alpha * X[i * xs0]; });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const int m,
                                           const ScalarType* KOKKOS_RESTRICT alpha, const int alphas0,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m),
                         [&](const int& i) { Y[i * ys0] += alpha[i * alphas0] * X[i * xs0]; });
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
        SerialAxpyInternal::invoke(n, alpha[i * alphas0], X + i * xs0, xs1, Y + i * ys0, ys1);
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int& j) {
        SerialAxpyInternal::invoke(m, alpha, alphas0, X + j * xs1, xs0, Y + j * ys1, ys0);
      });
    }
    // member.team_barrier();
    return 0;
  }
};

///
/// TeamVector Internal Impl
/// ========================
struct TeamVectorAxpyInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const int m, const ScalarType alpha,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m), [&](const int& i) { Y[i * ys0] += alpha * X[i * xs0]; });
    // member.team_barrier();
    return 0;
  }

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const int m,
                                           const ScalarType* KOKKOS_RESTRICT alpha, const int alphas0,
                                           const ValueType* KOKKOS_RESTRICT X, const int xs0,
                                           /* */ ValueType* KOKKOS_RESTRICT Y, const int ys0) {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m),
                         [&](const int& i) { Y[i * ys0] += alpha[i * alphas0] * X[i * xs0]; });
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
      Y[i * ys0 + j * ys1] += alpha[i * alphas0] * X[i * xs0 + j * xs1];
    });
    // member.team_barrier();
    return 0;
  }
};

///
/// Serial Impl
/// ===========

template <typename XViewType, typename YViewType, typename alphaViewType>
KOKKOS_INLINE_FUNCTION int SerialAxpy::invoke(const alphaViewType& alpha, const XViewType& X, const YViewType& Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert(Kokkos::is_view<XViewType>::value, "KokkosBatched::axpy: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<YViewType>::value, "KokkosBatched::axpy: YViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<alphaViewType>::value, "KokkosBatched::axpy: alphaViewType is not a Kokkos::View.");
  static_assert(XViewType::rank == 2, "KokkosBatched::axpy: XViewType must have rank 2.");
  static_assert(YViewType::rank == 2, "KokkosBatched::axpy: YViewType must have rank 2.");
  static_assert(alphaViewType::rank == 1, "KokkosBatched::axpy: alphaViewType must have rank 1.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
    Kokkos::printf(
        "KokkosBatched::axpy: Dimensions of X and Y do not match: X: %d x %d, "
        "Y: %d x %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
    return 1;
  }
  if (X.extent(0) != alpha.extent(0)) {
    Kokkos::printf(
        "KokkosBatched::axpy: First dimension of X and alpha do not match: X: "
        "%d x %d, alpha: %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)alpha.extent(0));
    return 1;
  }
#endif

  // No need to check if X.extent(0)==1 in the serial case as we don't
  // parallelize the kernel anyway.

  return SerialAxpyInternal::template invoke<typename alphaViewType::non_const_value_type,
                                             typename XViewType::non_const_value_type>(
      X.extent(0), X.extent(1), alpha.data(), alpha.stride_0(), X.data(), X.stride_0(), X.stride_1(), Y.data(),
      Y.stride_0(), Y.stride_1());
}

///
/// Team Impl
/// =========

template <typename MemberType>
template <typename XViewType, typename YViewType, typename alphaViewType>
KOKKOS_INLINE_FUNCTION int TeamAxpy<MemberType>::invoke(const MemberType& member, const alphaViewType& alpha,
                                                        const XViewType& X, const YViewType& Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert(Kokkos::is_view<XViewType>::value, "KokkosBatched::axpy: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<YViewType>::value, "KokkosBatched::axpy: YViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<alphaViewType>::value, "KokkosBatched::axpy: alphaViewType is not a Kokkos::View.");
  static_assert(XViewType::rank == 2, "KokkosBatched::axpy: XViewType must have rank 2.");
  static_assert(YViewType::rank == 2, "KokkosBatched::axpy: YViewType must have rank 2.");
  static_assert(alphaViewType::rank == 1, "KokkosBatched::axpy: alphaViewType must have rank 1.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
    Kokkos::printf(
        "KokkosBatched::axpy: Dimensions of X and Y do not match: X: %d x %d, "
        "Y: %d x %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
    return 1;
  }
  if (X.extent(0) != alpha.extent(0)) {
    Kokkos::printf(
        "KokkosBatched::axpy: First dimension of X and alpha do not match: X: "
        "%d x %d, alpha: %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)alpha.extent(0));
    return 1;
  }
#endif

  if (X.extent(0) == 1) {
    KokkosBlas::Experimental::axpy<MemberType>(member, alpha.data()[0], Kokkos::subview(X, 0, Kokkos::ALL),
                                               Kokkos::subview(Y, 0, Kokkos::ALL));
    return 0;
  }

  return TeamAxpyInternal::template invoke<MemberType, typename alphaViewType::non_const_value_type,
                                           typename XViewType::non_const_value_type>(
      member, X.extent(0), X.extent(1), alpha.data(), alpha.stride_0(), X.data(), X.stride_0(), X.stride_1(), Y.data(),
      Y.stride_0(), Y.stride_1());
}

///
/// TeamVector Impl
/// ===============

template <typename MemberType>
template <typename XViewType, typename YViewType, typename alphaViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorAxpy<MemberType>::invoke(const MemberType& member, const alphaViewType& alpha,
                                                              const XViewType& X, const YViewType& Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert(Kokkos::is_view<XViewType>::value, "KokkosBatched::axpy: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<YViewType>::value, "KokkosBatched::axpy: YViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<alphaViewType>::value, "KokkosBatched::axpy: alphaViewType is not a Kokkos::View.");
  static_assert(XViewType::rank == 2, "KokkosBatched::axpy: XViewType must have rank 2.");
  static_assert(YViewType::rank == 2, "KokkosBatched::axpy: YViewType must have rank 2.");
  static_assert(alphaViewType::rank == 1, "KokkosBatched::axpy: alphaViewType must have rank 1.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
    Kokkos::printf(
        "KokkosBatched::axpy: Dimensions of X and Y do not match: X: %d x %d, "
        "Y: %d x %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
    return 1;
  }
  if (X.extent(0) != alpha.extent(0)) {
    Kokkos::printf(
        "KokkosBatched::axpy: First dimension of X and alpha do not match: X: "
        "%d x %d, alpha: %d\n",
        (int)X.extent(0), (int)X.extent(1), (int)alpha.extent(0));
    return 1;
  }
#endif

  if (X.extent(0) == 1) {
    KokkosBlas::Experimental::axpy<MemberType>(member, alpha.data()[0], Kokkos::subview(X, 0, Kokkos::ALL),
                                               Kokkos::subview(Y, 0, Kokkos::ALL));
    return 0;
  }

  return TeamVectorAxpyInternal::invoke<MemberType, typename alphaViewType::non_const_value_type,
                                        typename XViewType::non_const_value_type, typename XViewType::array_layout>(
      member, X.extent(0), X.extent(1), alpha.data(), alpha.stride_0(), X.data(), X.stride_0(), X.stride_1(), Y.data(),
      Y.stride_0(), Y.stride_1());
}

}  // namespace KokkosBatched

#endif
