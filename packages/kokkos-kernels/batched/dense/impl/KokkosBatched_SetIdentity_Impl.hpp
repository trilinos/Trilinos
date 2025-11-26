// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SET_IDENTITY_IMPL_HPP
#define KOKKOSBATCHED_SET_IDENTITY_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_SetIdentity_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// ===========

template <typename AViewType>
KOKKOS_INLINE_FUNCTION int SerialSetIdentity::invoke(const AViewType &A) {
  return SerialSetIdentityInternal::invoke(A.extent(0), A.extent(1), A.data(), A.stride(0), A.stride(1));
}

///
/// Team Impl
/// =========

template <typename MemberType>
template <typename AViewType>
KOKKOS_INLINE_FUNCTION int TeamSetIdentity<MemberType>::invoke(const MemberType &member, const AViewType &A) {
  return TeamSetIdentityInternal::invoke(member, A.extent(0), A.extent(1), A.data(), A.stride(0), A.stride(1));
}

}  // end namespace KokkosBatched

#endif
