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
#ifndef __KOKKOSBATCHED_SET_IDENTITY_IMPL_HPP__
#define __KOKKOSBATCHED_SET_IDENTITY_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_SetIdentity_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// ===========

template <typename AViewType>
KOKKOS_INLINE_FUNCTION int SerialSetIdentity::invoke(const AViewType &A) {
  return SerialSetIdentityInternal::invoke(A.extent(0), A.extent(1), A.data(), A.stride_0(), A.stride_1());
}

///
/// Team Impl
/// =========

template <typename MemberType>
template <typename AViewType>
KOKKOS_INLINE_FUNCTION int TeamSetIdentity<MemberType>::invoke(const MemberType &member, const AViewType &A) {
  return TeamSetIdentityInternal::invoke(member, A.extent(0), A.extent(1), A.data(), A.stride_0(), A.stride_1());
}

}  // end namespace KokkosBatched

#endif
