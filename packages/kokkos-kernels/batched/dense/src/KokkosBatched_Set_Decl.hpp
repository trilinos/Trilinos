// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SET_DECL_HPP
#define KOKKOSBATCHED_SET_DECL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "impl/Kokkos_Error.hpp"

namespace KokkosBatched {
///
/// Serial Set
///

struct [[deprecated]] SerialSet {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A) {
    Kokkos::abort(
        "KokkosBatched::SerialSet is deprecated: use KokkosBlas::SerialSet "
        "instead");
    return 0;
  }  // namespace KokkosBatched
};

///
/// Team Set
///

template <typename MemberType>
struct [[deprecated]] TeamSet {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A) {
    Kokkos::abort(
        "KokkosBatched::TeamSet is deprecated: use KokkosBlas::TeamSet "
        "instead");
    return 0;
  }
};

///
/// TeamVector Set
///

template <typename MemberType>
struct [[deprecated]] TeamVectorSet {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A) {
    Kokkos::abort(
        "KokkosBatched::TeamVectorSet is deprecated: use "
        "KokkosBlas::TeamVectorSet instead");
    return 0;
  }
};

}  // namespace KokkosBatched

#endif
