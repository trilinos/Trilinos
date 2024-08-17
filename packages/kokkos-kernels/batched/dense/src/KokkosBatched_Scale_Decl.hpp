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
#ifndef __KOKKOSBATCHED_SCALE_DECL_HPP__
#define __KOKKOSBATCHED_SCALE_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "impl/Kokkos_Error.hpp"

namespace KokkosBatched {

///
/// Serial Scale
///

struct [[deprecated]] SerialScale {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A) {
    Kokkos::abort(
        "KokkosBatched::SerialScale is deprecated: use KokkosBlas::SerialScale "
        "instead");
    return 0;
  }  // namespace KokkosBatched
};

///
/// Team Scale
///

template <typename MemberType>
struct [[deprecated]] TeamScale {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A) {
    Kokkos::abort(
        "KokkosBatched::TeamScale is deprecated: use KokkosBlas::TeamScale "
        "instead");
    return 0;
  }
};

///
/// TeamVector Scale
///

template <typename MemberType>
struct [[deprecated]] TeamVectorScale {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A) {
    Kokkos::abort(
        "KokkosBatched::TeamVectorScale is deprecated: use "
        "KokkosBlas::TeamVectorScale instead");
    return 0;
  }
};

}  // namespace KokkosBatched

#endif
