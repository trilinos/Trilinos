// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_TRSV_TEAM_IMPL_HPP
#define KOKKOSBATCHED_TRSV_TEAM_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trsv_Team_Internal.hpp"

namespace KokkosBatched {
///
/// Team Impl
/// ===========

///
/// Implemented:
/// L/NT, U/NT, L/T, U/T
///
/// Not yet implemented
/// L/CT, U/CT

///
/// L/NT
///

template <typename MemberType, typename ArgDiag>
struct TeamTrsv<MemberType, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const bViewType &b) {
    return TeamTrsvInternalLower<Algo::Trsv::Unblocked>::invoke(
        member, ArgDiag::use_unit_diag, A.extent(0), alpha, A.data(), A.stride(0), A.stride(1), b.data(), b.stride(0));
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsv<MemberType, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const bViewType &b) {
    return TeamTrsvInternalLower<Algo::Trsv::Blocked>::invoke(
        member, ArgDiag::use_unit_diag, A.extent(0), alpha, A.data(), A.stride(0), A.stride(1), b.data(), b.stride(0));
  }
};

///
/// L/T
///

template <typename MemberType, typename ArgDiag>
struct TeamTrsv<MemberType, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const bViewType &b) {
    return TeamTrsvInternalUpper<Algo::Trsv::Unblocked>::invoke(
        member, ArgDiag::use_unit_diag, A.extent(1), alpha, A.data(), A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsv<MemberType, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType & /*member*/, const ScalarType alpha, const AViewType &A,
                                           const bViewType &b) {
    return TeamTrsvInternalUpper<Algo::Trsv::Blocked>::invoke(ArgDiag::use_unit_diag, A.extent(1), alpha, A.data(),
                                                              A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};

///
/// U/NT
///

template <typename MemberType, typename ArgDiag>
struct TeamTrsv<MemberType, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const bViewType &b) {
    return TeamTrsvInternalUpper<Algo::Trsv::Unblocked>::invoke(
        member, ArgDiag::use_unit_diag, A.extent(0), alpha, A.data(), A.stride(0), A.stride(1), b.data(), b.stride(0));
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsv<MemberType, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const bViewType &b) {
    return TeamTrsvInternalUpper<Algo::Trsv::Blocked>::invoke(
        member, ArgDiag::use_unit_diag, A.extent(0), alpha, A.data(), A.stride(0), A.stride(1), b.data(), b.stride(0));
  }
};

///
/// U/T
///

template <typename MemberType, typename ArgDiag>
struct TeamTrsv<MemberType, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const bViewType &b) {
    return TeamTrsvInternalLower<Algo::Trsv::Unblocked>::invoke(
        member, ArgDiag::use_unit_diag, A.extent(1), alpha, A.data(), A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsv<MemberType, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsv::Blocked> {
  template <typename ScalarType, typename AViewType, typename bViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const bViewType &b) {
    return TeamTrsvInternalLower<Algo::Trsv::Blocked>::invoke(
        member, ArgDiag::use_unit_diag, A.extent(1), alpha, A.data(), A.stride(1), A.stride(0), b.data(), b.stride(0));
  }
};
}  // namespace KokkosBatched

#endif
