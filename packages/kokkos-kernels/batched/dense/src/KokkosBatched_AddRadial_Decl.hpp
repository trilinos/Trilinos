// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_ADD_RADIAL_DECL_HPP
#define KOKKOSBATCHED_ADD_RADIAL_DECL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

///
/// This add tiny values on diagonals so the absolute values of diagonals become
/// larger
///

///
/// Serial AddRadial
///

struct SerialAddRadial {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType tiny, const AViewType &A);
};

///
/// Team Set
///

template <typename MemberType>
struct TeamAddRadial {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType tiny, const AViewType &A);
};

}  // namespace KokkosBatched

#include "KokkosBatched_AddRadial_Impl.hpp"

#endif
