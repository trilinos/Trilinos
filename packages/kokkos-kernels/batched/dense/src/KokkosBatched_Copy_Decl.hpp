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
#ifndef __KOKKOSBATCHED_COPY_DECL_HPP__
#define __KOKKOSBATCHED_COPY_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Copy
///

template <typename ArgTrans = Trans::NoTranspose, int rank = 2>
struct SerialCopy {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const BViewType &B);
};

///
/// Team Copy
///

template <typename MemberType, typename ArgTrans = Trans::NoTranspose, int rank = 2>
struct TeamCopy {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B);
};

///
/// TeamVector Copy
///

template <typename MemberType, typename ArgTrans = Trans::NoTranspose, int rank = 2>
struct TeamVectorCopy {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B);
};

///
/// Selective Interface
///
template <typename MemberType, typename ArgTrans, typename ArgMode, int rank = 2>
struct Copy {
  template <typename AViewType, typename BViewType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    int r_val = 0;
    if (std::is_same<ArgMode, Mode::Serial>::value) {
      r_val = SerialCopy<ArgTrans, rank>::invoke(A, B);
    } else if (std::is_same<ArgMode, Mode::Team>::value) {
      r_val = TeamCopy<MemberType, ArgTrans, rank>::invoke(member, A, B);
    } else if (std::is_same<ArgMode, Mode::TeamVector>::value) {
      r_val = TeamVectorCopy<MemberType, ArgTrans, rank>::invoke(member, A, B);
    }
    return r_val;
  }
};

}  // namespace KokkosBatched

#include "KokkosBatched_Copy_Impl.hpp"

#define KOKKOSBATCHED_SERIAL_COPY_MATRIX_NO_TRANSPOSE_INTERNAL_INVOKE(M, N, A, AS0, AS1, B, BS0, BS1) \
  KokkosBatched::SerialCopyInternal ::invoke(M, N, A, AS0, AS1, B, BS0, BS1)

#define KOKKOSBATCHED_TEAM_COPY_MATRIX_NO_TRANSPOSE_INTERNAL_INVOKE(MEMBER, M, N, A, AS0, AS1, B, BS0, BS1) \
  KokkosBatched::TeamCopyInternal ::invoke(MEMBER, M, N, A, AS0, AS1, B, BS0, BS1)

#define KOKKOSBATCHED_SERIAL_COPY_VECTOR_INTERNAL_INVOKE(M, A, AS, B, BS) \
  KokkosBatched::SerialCopyInternal ::invoke(M, A, AS, B, BS)

#define KOKKOSBATCHED_TEAM_COPY_VECTOR_NO_TRANSPOSE_INTERNAL_INVOKE(MEMBER, M, A, AS, B, BS) \
  KokkosBatched::TeamCopyInternal ::invoke(MEMBER, M, A, AS, B, BS)

#define KOKKOSBATCHED_COPY_VECTOR_NO_TRANSPOSE_INTERNAL_INVOKE(MODETYPE, MEMBER, M, A, AS, B, BS) \
  if (std::is_same<MODETYPE, KokkosBatched::Mode::Serial>::value) {                               \
    KOKKOSBATCHED_SERIAL_COPY_VECTOR_INTERNAL_INVOKE(M, A, AS, B, BS);                            \
  } else if (std::is_same<MODETYPE, KokkosBatched::Mode::Team>::value) {                          \
    KOKKOSBATCHED_TEAM_COPY_VECTOR_NO_TRANSPOSE_INTERNAL_INVOKE(MEMBER, M, A, AS, B, BS);         \
  }

#endif
