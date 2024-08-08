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
#ifndef __KOKKOSBATCHED_DOT_INTERNAL_HPP__
#define __KOKKOSBATCHED_DOT_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBlas1_team_dot.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================

struct SerialDotInternal {
  // i \in [0,m)
  // C = conj(A(:))*B(:)
  template <typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const int m, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C) {
    using ats = Kokkos::ArithTraits<ValueType>;
    C[0]      = ValueType(0);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) {
      const int idx_a = i * as0, idx_b = i * bs0;
      C[0] += ats::conj(A[idx_a]) * B[idx_b];
    }
    return 0;
  }

  // j \in [0,n), i \in [0,m)
  // C(j) = conj(A(:,j))*B(:,j)
  template <typename ValueType, typename MagnitudeType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int n, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           const int as1, const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                           const int bs1,
                                           /* */ MagnitudeType *KOKKOS_RESTRICT C, const int cs) {
    for (int j = 0; j < n; ++j) invoke(m, A + j * as1, as0, B + j * bs1, bs0, C + j * cs);
    return 0;
  }
};

///
/// Team Internal Impl
/// ========================

// i \in [0,m)
// C = conj(A(:))*B(:)
struct TeamDotInternal {
  template <typename MemberType, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C) {
    using ats = Kokkos::ArithTraits<ValueType>;
    ValueType t(0);
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(member, m),
        [&](const int &i, ValueType &update) {
          const int idx_a = i * as0, idx_b = i * bs0;
          update += ats::conj(A[idx_a]) * B[idx_b];
        },
        t);
    Kokkos::single(Kokkos::PerThread(member), [&]() { C[0] = t; });
    return 0;
  }

  // j \in [0,n), i \in [0,m)
  // C(j) = conj(A(:,j))*B(:,j)
  template <typename MemberType, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C, const int cs) {
    using ats = Kokkos::ArithTraits<ValueType>;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &j) {
      ValueType t(0);
      const ValueType *KOKKOS_RESTRICT A_at_j = A + j * as1;
      const ValueType *KOKKOS_RESTRICT B_at_j = B + j * bs1;
      for (int i = 0; i < m; ++i) {
        const int idx_a = i * as0, idx_b = i * bs0;
        t += ats::conj(A_at_j[idx_a]) * B_at_j[idx_b];
      }
      Kokkos::single(Kokkos::PerThread(member), [&]() { C[j * cs] = t; });
    });
    return 0;
  }
};

///
/// TeamVector Internal Impl
/// ========================

// i \in [0,m)
// C = conj(A(:))*B(:)
struct TeamVectorDotInternal {
  template <typename MemberType, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C) {
    using ats = Kokkos::ArithTraits<ValueType>;
    ValueType t(0);
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(member, m),
        [&](const int &i, ValueType &update) {
          const int idx_a = i * as0, idx_b = i * bs0;
          update += ats::conj(A[idx_a]) * B[idx_b];
        },
        t);
    Kokkos::single(Kokkos::PerThread(member), [&]() { C[0] = t; });
    return 0;
  }

  // j \in [0,n), i \in [0,m)
  // C(j) = conj(A(:,j))*B(:,j)
  template <typename MemberType, typename ValueType, typename MagnitudeType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const int m, const int n,
                                                const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                                const ValueType *KOKKOS_RESTRICT B, const int bs0, const int bs1,
                                                /* */ MagnitudeType *KOKKOS_RESTRICT C, const int cs) {
    using ats = Kokkos::ArithTraits<ValueType>;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const int &j) {
      ValueType t(0);
      const ValueType *KOKKOS_RESTRICT A_at_j = A + j * as1;
      const ValueType *KOKKOS_RESTRICT B_at_j = B + j * bs1;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(member, m),
          [&](const int &i, ValueType &update) {
            const int idx_a = i * as0, idx_b = i * bs0;
            update += ats::conj(A_at_j[idx_a]) * B_at_j[idx_b];
          },
          t);
      Kokkos::single(Kokkos::PerThread(member), [&]() { C[j * cs] = t; });
    });
    return 0;
  }
};

///
/// Serial Impl
/// ===========

template <>
struct SerialDot<Trans::Transpose> {
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &X, const YViewType &Y, const NormViewType &dot) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<XViewType>::value, "KokkosBatched::dot: XViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YViewType>::value, "KokkosBatched::dot: YViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<NormViewType>::value, "KokkosBatched::dot: NormViewType is not a Kokkos::View.");
    static_assert(XViewType::rank == 2, "KokkosBatched::dot: XViewType must have rank 2.");
    static_assert(YViewType::rank == 2, "KokkosBatched::dot: YViewType must have rank 2.");
    static_assert(NormViewType::rank == 1, "KokkosBatched::dot: NormViewType must have rank 1.");

    // Check compatibility of dimensions at run time.
    if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::dot: Dimensions of X and Y do not match: X: %d x %d, "
          "Y: %d x %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
      return 1;
    }
    if (X.extent(1) != dot.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::dot: Second dimension of X and alpha do not match: "
          "X: "
          "%d x %d, dot: %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)dot.extent(0));
      return 1;
    }
#endif
    return SerialDotInternal::template invoke<typename XViewType::non_const_value_type,
                                              typename NormViewType::non_const_value_type>(
        X.extent(0), X.extent(1), X.data(), X.stride_0(), X.stride_1(), Y.data(), Y.stride_0(), Y.stride_1(),
        dot.data(), dot.stride_0());
  }
};

template <>
struct SerialDot<Trans::NoTranspose> {
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &X, const YViewType &Y, const NormViewType &dot) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<XViewType>::value, "KokkosBatched::dot: XViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YViewType>::value, "KokkosBatched::dot: YViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<NormViewType>::value, "KokkosBatched::dot: NormViewType is not a Kokkos::View.");
    static_assert(XViewType::rank == 2, "KokkosBatched::dot: XViewType must have rank 2.");
    static_assert(YViewType::rank == 2, "KokkosBatched::dot: YViewType must have rank 2.");
    static_assert(NormViewType::rank == 1, "KokkosBatched::dot: NormViewType must have rank 1.");

    // Check compatibility of dimensions at run time.
    if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::dot: Dimensions of X and Y do not match: X: %d x %d, "
          "Y: %d x %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
      return 1;
    }
    if (X.extent(0) != dot.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::dot: First dimension of X and alpha do not match: X: "
          "%d x %d, dot: %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)dot.extent(0));
      return 1;
    }
#endif
    return SerialDotInternal::template invoke<typename XViewType::non_const_value_type,
                                              typename NormViewType::non_const_value_type>(
        X.extent(1), X.extent(0), X.data(), X.stride_1(), X.stride_0(), Y.data(), Y.stride_1(), Y.stride_0(),
        dot.data(), dot.stride_0());
  }
};

///
/// Team Impl
/// ===============

template <typename MemberType>
struct TeamDot<MemberType, Trans::Transpose> {
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const YViewType &Y,
                                           const NormViewType &dot) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<XViewType>::value, "KokkosBatched::dot: XViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YViewType>::value, "KokkosBatched::dot: YViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<NormViewType>::value, "KokkosBatched::dot: NormViewType is not a Kokkos::View.");
    static_assert(XViewType::rank == 2, "KokkosBatched::dot: XViewType must have rank 2.");
    static_assert(YViewType::rank == 2, "KokkosBatched::dot: YViewType must have rank 2.");
    static_assert(NormViewType::rank == 1, "KokkosBatched::dot: NormViewType must have rank 1.");

    // Check compatibility of dimensions at run time.
    if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::dot: Dimensions of X and Y do not match: X: %d x %d, "
          "Y: %d x %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
      return 1;
    }
    if (X.extent(1) != dot.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::dot: Second dimension of X and alpha do not match: "
          "X: "
          "%d x %d, dot: %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)dot.extent(0));
      return 1;
    }
#endif

    if (X.extent(1) == 1) {
      dot(0) =
          KokkosBlas::Experimental::dot(member, Kokkos::subview(X, Kokkos::ALL, 0), Kokkos::subview(Y, Kokkos::ALL, 0));
      return 0;
    }

    return TeamDotInternal::template invoke<MemberType, typename XViewType::non_const_value_type,
                                            typename NormViewType::non_const_value_type>(
        member, X.extent(0), X.extent(1), X.data(), X.stride_0(), X.stride_1(), Y.data(), Y.stride_0(), Y.stride_1(),
        dot.data(), dot.stride_0());
  }
};

template <typename MemberType>
struct TeamDot<MemberType, Trans::NoTranspose> {
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const YViewType &Y,
                                           const NormViewType &dot) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<XViewType>::value, "KokkosBatched::dot: XViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YViewType>::value, "KokkosBatched::dot: YViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<NormViewType>::value, "KokkosBatched::dot: NormViewType is not a Kokkos::View.");
    static_assert(XViewType::rank == 2, "KokkosBatched::dot: XViewType must have rank 2.");
    static_assert(YViewType::rank == 2, "KokkosBatched::dot: YViewType must have rank 2.");
    static_assert(NormViewType::rank == 1, "KokkosBatched::dot: NormViewType must have rank 1.");

    // Check compatibility of dimensions at run time.
    if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::dot: Dimensions of X and Y do not match: X: %d x %d, "
          "Y: %d x %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
      return 1;
    }
    if (X.extent(0) != dot.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::dot: First dimension of X and alpha do not match: X: "
          "%d x %d, dot: %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)dot.extent(0));
      return 1;
    }
#endif

    if (X.extent(0) == 1) {
      dot(0) =
          KokkosBlas::Experimental::dot(member, Kokkos::subview(X, 0, Kokkos::ALL), Kokkos::subview(Y, 0, Kokkos::ALL));
      return 0;
    }

    return TeamDotInternal::template invoke<MemberType, typename XViewType::non_const_value_type,
                                            typename NormViewType::non_const_value_type>(
        member, X.extent(1), X.extent(0), X.data(), X.stride_1(), X.stride_0(), Y.data(), Y.stride_1(), Y.stride_0(),
        dot.data(), dot.stride_0());
  }
};

///
/// TeamVector Impl
/// ===============

template <typename MemberType>
struct TeamVectorDot<MemberType, Trans::Transpose> {
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const YViewType &Y,
                                           const NormViewType &dot) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<XViewType>::value, "KokkosBatched::dot: XViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YViewType>::value, "KokkosBatched::dot: YViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<NormViewType>::value, "KokkosBatched::dot: NormViewType is not a Kokkos::View.");
    static_assert(XViewType::rank == 2, "KokkosBatched::dot: XViewType must have rank 2.");
    static_assert(YViewType::rank == 2, "KokkosBatched::dot: YViewType must have rank 2.");
    static_assert(NormViewType::rank == 1, "KokkosBatched::dot: NormViewType must have rank 1.");

    // Check compatibility of dimensions at run time.
    if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::dot: Dimensions of X and Y do not match: X: %d x %d, "
          "Y: %d x %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
      return 1;
    }
    if (X.extent(1) != dot.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::dot: Second dimension of X and alpha do not match: "
          "X: "
          "%d x %d, dot: %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)dot.extent(0));
      return 1;
    }
#endif

    if (X.extent(1) == 1) {
      dot(0) =
          KokkosBlas::Experimental::dot(member, Kokkos::subview(X, Kokkos::ALL, 0), Kokkos::subview(Y, Kokkos::ALL, 0));
      return 0;
    }

    return TeamVectorDotInternal::template invoke<MemberType, typename XViewType::non_const_value_type,
                                                  typename NormViewType::non_const_value_type>(
        member, X.extent(0), X.extent(1), X.data(), X.stride_0(), X.stride_1(), Y.data(), Y.stride_0(), Y.stride_1(),
        dot.data(), dot.stride_0());
  }
};

template <typename MemberType>
struct TeamVectorDot<MemberType, Trans::NoTranspose> {
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const YViewType &Y,
                                           const NormViewType &dot) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<XViewType>::value, "KokkosBatched::dot: XViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<YViewType>::value, "KokkosBatched::dot: YViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<NormViewType>::value, "KokkosBatched::dot: NormViewType is not a Kokkos::View.");
    static_assert(XViewType::rank == 2, "KokkosBatched::dot: XViewType must have rank 2.");
    static_assert(YViewType::rank == 2, "KokkosBatched::dot: YViewType must have rank 2.");
    static_assert(NormViewType::rank == 1, "KokkosBatched::dot: NormViewType must have rank 1.");

    // Check compatibility of dimensions at run time.
    if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::dot: Dimensions of X and Y do not match: X: %d x %d, "
          "Y: %d x %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
      return 1;
    }
    if (X.extent(0) != dot.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::dot: First dimension of X and alpha do not match: X: "
          "%d x %d, dot: %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)dot.extent(0));
      return 1;
    }
#endif

    if (X.extent(0) == 1) {
      dot(0) =
          KokkosBlas::Experimental::dot(member, Kokkos::subview(X, 0, Kokkos::ALL), Kokkos::subview(Y, 0, Kokkos::ALL));
      return 0;
    }

    return TeamVectorDotInternal::template invoke<MemberType, typename XViewType::non_const_value_type,
                                                  typename NormViewType::non_const_value_type>(
        member, X.extent(1), X.extent(0), X.data(), X.stride_1(), X.stride_0(), Y.data(), Y.stride_1(), Y.stride_0(),
        dot.data(), dot.stride_0());
  }
};

}  // end namespace KokkosBatched

#endif
