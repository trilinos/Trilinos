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
#ifndef __KOKKOSBATCHED_TRSM_TEAM_IMPL_HPP__
#define __KOKKOSBATCHED_TRSM_TEAM_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trsm_Team_Internal.hpp"

namespace KokkosBatched {

///
/// Team Impl
/// =========

///
/// L/L/NT
///
/// B := inv(tril(A)) (alpha*B)
/// A(m x m), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                    B.extent(1), alpha, A.data(), A.stride_0(),
                                                                    A.stride_1(), B.data(), B.stride_0(), B.stride_1());
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                  B.extent(1), alpha, A.data(), A.stride_0(),
                                                                  A.stride_1(), B.data(), B.stride_0(), B.stride_1());
  }
};

///
/// R/U/NT
///
/// B := (alpha*B) inv(triu(A))
/// A(n x n), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(1),
                                                                    B.extent(0), alpha, A.data(), A.stride_1(),
                                                                    A.stride_0(), B.data(), B.stride_1(), B.stride_0());
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(1),
                                                                  B.extent(0), alpha, A.data(), A.stride_1(),
                                                                  A.stride_0(), B.data(), B.stride_1(), B.stride_0());
  }
};

///
/// R/L/NT
///
/// B := (alpha*B) inv(tril(A))
/// A(n x n), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(1),
                                                                    B.extent(0), alpha, A.data(), A.stride_1(),
                                                                    A.stride_0(), B.data(), B.stride_1(), B.stride_0());
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(1),
                                                                  B.extent(0), alpha, A.data(), A.stride_1(),
                                                                  A.stride_0(), B.data(), B.stride_1(), B.stride_0());
  }
};

///
/// R/U/T
///
/// B := (alpha*B) inv(triu(A))
/// A(n x n), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(1),
                                                                    B.extent(0), alpha, A.data(), A.stride_0(),
                                                                    A.stride_1(), B.data(), B.stride_1(), B.stride_0());
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(1),
                                                                  B.extent(0), alpha, A.data(), A.stride_0(),
                                                                  A.stride_1(), B.data(), B.stride_1(), B.stride_0());
  }
};

///
/// L/U/NT
///
/// B := inv(triu(A)) (alpha*B)
/// A(m x m), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                    B.extent(1), alpha, A.data(), A.stride_0(),
                                                                    A.stride_1(), B.data(), B.stride_0(), B.stride_1());
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                  B.extent(1), alpha, A.data(), A.stride_0(),
                                                                  A.stride_1(), B.data(), B.stride_0(), B.stride_1());
  }
};

///
/// L/L/T
///
/// B := inv(tril(AT)) (alpha*B)
/// A(m x m), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                    B.extent(1), alpha, A.data(), A.stride_1(),
                                                                    A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                  B.extent(1), alpha, A.data(), A.stride_1(),
                                                                  A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};

///
/// L/U/T
///
/// B := inv(triu(AT)) (alpha*B)
/// A(m x m), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                    B.extent(1), alpha, A.data(), A.stride_1(),
                                                                    A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    return TeamTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                  B.extent(1), alpha, A.data(), A.stride_1(),
                                                                  A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};

}  // namespace KokkosBatched

#endif
