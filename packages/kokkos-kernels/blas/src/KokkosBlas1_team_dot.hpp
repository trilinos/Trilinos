// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
