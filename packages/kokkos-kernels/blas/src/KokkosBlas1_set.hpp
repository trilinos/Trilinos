// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_SET_HPP_
#define KOKKOSBLAS1_SET_HPP_

#include <KokkosBlas1_set_impl.hpp>

namespace KokkosBlas {

///
/// Serial Set
///

struct SerialSet {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A) {
    if constexpr (AViewType::rank() == 1)
      return Impl::SerialSetInternal::invoke(A.extent(0), alpha, A.data(), A.stride(0));
    else
      return Impl::SerialSetInternal::invoke(A.extent(0), A.extent(1), alpha, A.data(), A.stride(0), A.stride(1));
  }
};

///
/// Team Set
///

template <typename MemberType>
struct TeamSet {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A) {
    if constexpr (AViewType::rank() == 1)
      return Impl::TeamSetInternal::invoke(member, A.extent(0), alpha, A.data(), A.stride(0));
    else
      return Impl::TeamSetInternal::invoke(member, A.extent(0), A.extent(1), alpha, A.data(), A.stride(0), A.stride(1));
  }
};

///
/// TeamVector Set
///

template <typename MemberType>
struct TeamVectorSet {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A) {
    if constexpr (AViewType::rank() == 1)
      return Impl::TeamVectorSetInternal::invoke(member, A.extent(0), alpha, A.data(), A.stride(0));
    else
      return Impl::TeamVectorSetInternal::invoke(member, A.extent(0), A.extent(1), alpha, A.data(), A.stride(0),
                                                 A.stride(1));
  }
};

}  // namespace KokkosBlas

#endif
