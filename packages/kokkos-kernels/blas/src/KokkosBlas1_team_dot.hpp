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

#ifndef KOKKOSBLAS1_TEAM_DOT_HPP_
#define KOKKOSBLAS1_TEAM_DOT_HPP_

#include <KokkosBlas1_team_dot_spec.hpp>

namespace KokkosBlas {
namespace Experimental {

template <class TeamType, class XVector, class YVector>
typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type
    KOKKOS_INLINE_FUNCTION
    dot(const TeamType& team, const XVector& x, const YVector& y) {
  return Impl::TeamDot<TeamType, XVector, YVector>::team_dot(team, x, y);
}

}  // namespace Experimental
}  // namespace KokkosBlas

#endif
