// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_QR_TEAMVECTOR_IMPL_HPP
#define KOKKOSBATCHED_QR_TEAMVECTOR_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_QR_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Impl
/// ===============

template <typename MemberType>
struct TeamVectorQR<MemberType, Algo::QR::Unblocked> {
  template <typename AViewType, typename tViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const wViewType &w) {
    return TeamVectorQR_Internal::invoke(member, A.extent(0), A.extent(1), A.data(), A.stride(0), A.stride(1), t.data(),
                                         t.stride(0), w.data());
  }
};

}  // namespace KokkosBatched

#endif
