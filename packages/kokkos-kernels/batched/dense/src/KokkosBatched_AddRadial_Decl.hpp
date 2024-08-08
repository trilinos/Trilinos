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
#ifndef __KOKKOSBATCHED_ADD_RADIAL_DECL_HPP__
#define __KOKKOSBATCHED_ADD_RADIAL_DECL_HPP__

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

#endif
