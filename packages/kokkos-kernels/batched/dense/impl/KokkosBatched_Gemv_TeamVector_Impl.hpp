// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_GEMV_TEAMVECTOR_IMPL_HPP
#define KOKKOSBATCHED_GEMV_TEAMVECTOR_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Gemv_TeamVector_Internal.hpp"
#include "KokkosBlas2_team_gemv.hpp"

namespace KokkosBatched {

///
/// Team Impl
/// =========

///
/// Implemented:
/// NT, T
///
/// Not yet implemented
/// CT

///
/// NT
///

template <typename MemberType>
struct TeamVectorGemv<MemberType, Trans::NoTranspose, Algo::Gemv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const xViewType &x, const ScalarType beta, const yViewType &y) {
    static_assert(AViewType::rank == 3,
                  "Batched TeamVectorGemv requires rank-3 A matrix (use "
                  "KokkosBlas::TeamVectorGemv for regular rank-2 matrix)");
    if (A.extent(0) == 1) {
      KokkosBlas::TeamVectorGemv<MemberType, Trans::NoTranspose, Algo::Gemv::Unblocked>::invoke(
          member, alpha, Kokkos::subview(A, 0, Kokkos::ALL, Kokkos::ALL), Kokkos::subview(x, 0, Kokkos::ALL), beta,
          Kokkos::subview(y, 0, Kokkos::ALL));
      return 0;
    }
    return TeamVectorGemvInternal<Algo::Gemv::Unblocked>::template invoke<
        MemberType, ScalarType, typename AViewType::array_layout, typename AViewType::non_const_value_type>(
        member, A.extent(0), A.extent(1), A.extent(2), alpha, A.data(), A.stride(0), A.stride(1), A.stride(2), x.data(),
        x.stride(0), x.stride(1), beta, y.data(), y.stride(0), y.stride(1));
  }
};

template <typename MemberType>
struct TeamVectorGemv<MemberType, Trans::NoTranspose, Algo::Gemv::Blocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType & /*member*/, const ScalarType /*alpha*/,
                                           const AViewType & /*A*/, const xViewType & /*x*/, const ScalarType /*beta*/,
                                           const yViewType & /*y*/) {
    static_assert(AViewType::rank == 3,
                  "Batched TeamVectorGemv requires rank-3 A matrix (use "
                  "KokkosBlas::TeamVectorGemv for regular rank-2 matrix)");
    Kokkos::abort(
        "KokkosBatched::TeamVectorGemv<Algo::Gemv::Blocked> for rank-3 matrix "
        "is NOT implemented");
  }
};

///
/// T
///

template <typename MemberType>
struct TeamVectorGemv<MemberType, Trans::Transpose, Algo::Gemv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const xViewType &x, const ScalarType beta, const yViewType &y) {
    static_assert(AViewType::rank == 3,
                  "Batched TeamVectorGemv requires rank-3 A matrix (use "
                  "KokkosBlas::TeamVectorGemv for regular rank-2 matrix)");
    return TeamVectorGemvInternal<Algo::Gemv::Unblocked>::template invoke<
        MemberType, ScalarType, typename AViewType::array_layout, typename AViewType::non_const_value_type>(
        member, A.extent(0), A.extent(2), A.extent(1), alpha, A.data(), A.stride(0), A.stride(2), A.stride(1), x.data(),
        x.stride(0), x.stride(1), beta, y.data(), y.stride(0), y.stride(1));
  }
};

template <typename MemberType>
struct TeamVectorGemv<MemberType, Trans::Transpose, Algo::Gemv::Blocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType & /*member*/, const ScalarType /*alpha*/,
                                           const AViewType & /*A*/, const xViewType & /*x*/, const ScalarType /*beta*/,
                                           const yViewType & /*y*/) {
    static_assert(AViewType::rank == 3,
                  "Batched TeamVectorGemv requires rank-3 A matrix (use "
                  "KokkosBlas::TeamVectorGemv for regular rank-2 matrix)");
    Kokkos::abort(
        "KokkosBatched::TeamVectorGemv<Algo::Gemv::Blocked> for rank-3 matrix "
        "is NOT implemented");
  }
};

}  // namespace KokkosBatched

#endif
