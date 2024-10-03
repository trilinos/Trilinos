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
#ifndef __KOKKOSBATCHED_COPY_IMPL_HPP__
#define __KOKKOSBATCHED_COPY_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Copy_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// ===========

template <>
template <typename AViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialCopy<Trans::NoTranspose, 1>::invoke(const AViewType &A, const BViewType &B) {
  return SerialCopyInternal::invoke(A.extent(0), A.data(), A.stride_0(), B.data(), B.stride_0());
}

template <>
template <typename AViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialCopy<Trans::Transpose, 1>::invoke(const AViewType &A, const BViewType &B) {
  return SerialCopyInternal::invoke(A.extent(0), A.data(), A.stride_0(), B.data(), B.stride_0());
}

template <>
template <typename AViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialCopy<Trans::NoTranspose, 2>::invoke(const AViewType &A, const BViewType &B) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert(Kokkos::is_view<AViewType>::value, "KokkosBatched::copy: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<BViewType>::value, "KokkosBatched::copy: BViewType is not a Kokkos::View.");
  static_assert(AViewType::rank == 2, "KokkosBatched::copy: AViewType must have rank 2.");
  static_assert(BViewType::rank == 2, "KokkosBatched::copy: BViewType must have rank 2.");

  // Check compatibility of dimensions at run time.
  if (A.extent(0) != B.extent(0) || A.extent(1) != B.extent(1)) {
    Kokkos::printf(
        "KokkosBatched::copy: Dimensions of A and B do not match: A: %d x %d, "
        "B: %d x %d\n",
        (int)A.extent(0), (int)A.extent(1), (int)B.extent(0), (int)B.extent(1));
    return 1;
  }
#endif
  return SerialCopyInternal::invoke(A.extent(0), A.extent(1), A.data(), A.stride_0(), A.stride_1(), B.data(),
                                    B.stride_0(), B.stride_1());
}

template <>
template <typename AViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialCopy<Trans::Transpose, 2>::invoke(const AViewType &A, const BViewType &B) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert(Kokkos::is_view<AViewType>::value, "KokkosBatched::copy: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<BViewType>::value, "KokkosBatched::copy: BViewType is not a Kokkos::View.");
  static_assert(AViewType::rank == 2, "KokkosBatched::copy: AViewType must have rank 2.");
  static_assert(BViewType::rank == 2, "KokkosBatched::copy: BViewType must have rank 2.");

  // Check compatibility of dimensions at run time.
  if (A.extent(0) != B.extent(0) || A.extent(1) != B.extent(1)) {
    Kokkos::printf(
        "KokkosBatched::copy: Dimensions of A and B do not match: A: %d x %d, "
        "B: %d x %d\n",
        (int)A.extent(0), (int)A.extent(1), (int)B.extent(0), (int)B.extent(1));
    return 1;
  }
#endif
  return SerialCopyInternal::invoke(A.extent(1), A.extent(0), A.data(), A.stride_1(), A.stride_0(), B.data(),
                                    B.stride_0(), B.stride_1());
}

///
/// Team Impl
/// =========

template <typename MemberType>
struct TeamCopy<MemberType, Trans::NoTranspose, 1> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    return TeamCopyInternal::invoke(member, A.extent(0), A.data(), A.stride_0(), B.data(), B.stride_0());
  }
};

template <typename MemberType>
struct TeamCopy<MemberType, Trans::Transpose, 1> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    return TeamCopyInternal::invoke(member, A.extent(0), A.data(), A.stride_0(), B.data(), B.stride_0());
  }
};

template <typename MemberType>
struct TeamCopy<MemberType, Trans::NoTranspose, 2> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<AViewType>::value, "KokkosBatched::copy: AViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<BViewType>::value, "KokkosBatched::copy: BViewType is not a Kokkos::View.");
    static_assert(AViewType::rank == 2, "KokkosBatched::copy: AViewType must have rank 2.");
    static_assert(BViewType::rank == 2, "KokkosBatched::copy: BViewType must have rank 2.");

    // Check compatibility of dimensions at run time.
    if (A.extent(0) != B.extent(0) || A.extent(1) != B.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::copy: Dimensions of A and B do not match: A: %d x "
          "%d, "
          "B: %d x %d\n",
          (int)A.extent(0), (int)A.extent(1), (int)B.extent(0), (int)B.extent(1));
      return 1;
    }
#endif
    if (A.extent(0) == 1) {
      return TeamCopy<MemberType, Trans::NoTranspose, 1>::invoke(member, Kokkos::subview(A, 0, Kokkos::ALL),
                                                                 Kokkos::subview(B, 0, Kokkos::ALL));
    }
    return TeamCopyInternal::invoke(member, A.extent(0), A.extent(1), A.data(), A.stride_0(), A.stride_1(), B.data(),
                                    B.stride_0(), B.stride_1());
  }
};

template <typename MemberType>
struct TeamCopy<MemberType, Trans::Transpose, 2> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<AViewType>::value, "KokkosBatched::copy: AViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<BViewType>::value, "KokkosBatched::copy: BViewType is not a Kokkos::View.");
    static_assert(AViewType::rank == 2, "KokkosBatched::copy: AViewType must have rank 2.");
    static_assert(BViewType::rank == 2, "KokkosBatched::copy: BViewType must have rank 2.");

    // Check compatibility of dimensions at run time.
    if (A.extent(0) != B.extent(0) || A.extent(1) != B.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::copy: Dimensions of A and B do not match: A: %d x "
          "%d, "
          "B: %d x %d\n",
          (int)A.extent(0), (int)A.extent(1), (int)B.extent(0), (int)B.extent(1));
      return 1;
    }
#endif
    if (A.extent(1) == 1) {
      return TeamCopy<MemberType, Trans::Transpose, 1>::invoke(member, Kokkos::subview(A, Kokkos::ALL, 0),
                                                               Kokkos::subview(B, Kokkos::ALL, 0));
    }
    return TeamCopyInternal::invoke(member, A.extent(1), A.extent(0), A.data(), A.stride_1(), A.stride_0(), B.data(),
                                    B.stride_0(), B.stride_1());
  }
};

///
/// TeamVector Impl
/// =========

template <typename MemberType>
struct TeamVectorCopy<MemberType, Trans::NoTranspose, 1> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    return TeamVectorCopyInternal::invoke(member, A.extent(0), A.data(), A.stride_0(), B.data(), B.stride_0());
  }
};

template <typename MemberType>
struct TeamVectorCopy<MemberType, Trans::Transpose, 1> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    return TeamVectorCopyInternal::invoke(member, A.extent(0), A.data(), A.stride_0(), B.data(), B.stride_0());
  }
};

template <typename MemberType>
struct TeamVectorCopy<MemberType, Trans::NoTranspose, 2> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<AViewType>::value, "KokkosBatched::copy: AViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<BViewType>::value, "KokkosBatched::copy: BViewType is not a Kokkos::View.");
    static_assert(AViewType::rank == 2, "KokkosBatched::copy: AViewType must have rank 2.");
    static_assert(BViewType::rank == 2, "KokkosBatched::copy: BViewType must have rank 2.");

    // Check compatibility of dimensions at run time.
    if (A.extent(0) != B.extent(0) || A.extent(1) != B.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::copy: Dimensions of A and B do not match: A: %d x "
          "%d, "
          "B: %d x %d\n",
          (int)A.extent(0), (int)A.extent(1), (int)B.extent(0), (int)B.extent(1));
      return 1;
    }
#endif
    if (A.extent(0) == 1) {
      return TeamVectorCopy<MemberType, Trans::NoTranspose, 1>::invoke(member, Kokkos::subview(A, 0, Kokkos::ALL),
                                                                       Kokkos::subview(B, 0, Kokkos::ALL));
    }
    return TeamVectorCopyInternal::invoke(member, A.extent(0), A.extent(1), A.data(), A.stride_0(), A.stride_1(),
                                          B.data(), B.stride_0(), B.stride_1());
  }
};

template <typename MemberType>
struct TeamVectorCopy<MemberType, Trans::Transpose, 2> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<AViewType>::value, "KokkosBatched::copy: AViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<BViewType>::value, "KokkosBatched::copy: BViewType is not a Kokkos::View.");
    static_assert(AViewType::rank == 2, "KokkosBatched::copy: AViewType must have rank 2.");
    static_assert(BViewType::rank == 2, "KokkosBatched::copy: BViewType must have rank 2.");

    // Check compatibility of dimensions at run time.
    if (A.extent(0) != B.extent(0) || A.extent(1) != B.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::copy: Dimensions of A and B do not match: A: %d x "
          "%d, "
          "B: %d x %d\n",
          (int)A.extent(0), (int)A.extent(1), (int)B.extent(0), (int)B.extent(1));
      return 1;
    }
#endif
    if (A.extent(1) == 1) {
      return TeamVectorCopy<MemberType, Trans::NoTranspose, 1>::invoke(member, Kokkos::subview(A, Kokkos::ALL, 0),
                                                                       Kokkos::subview(B, Kokkos::ALL, 0));
    }
    return TeamVectorCopyInternal::invoke(member, A.extent(1), A.extent(0), A.data(), A.stride_1(), A.stride_0(),
                                          B.data(), B.stride_0(), B.stride_1());
  }
};

}  // end namespace KokkosBatched

#endif
