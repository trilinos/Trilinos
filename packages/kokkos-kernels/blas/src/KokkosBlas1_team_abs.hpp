// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_TEAM_ABS_HPP_
#define KOKKOSBLAS1_TEAM_ABS_HPP_

#include <KokkosBlas1_team_abs_spec.hpp>

namespace KokkosBlas {
namespace Experimental {

template <class TeamType, class RVector, class XVector>
void KOKKOS_INLINE_FUNCTION abs(const TeamType& team, const RVector& r, const XVector& x) {
  Impl::TeamAbs<TeamType, RVector, XVector>::team_abs(team, r, x);
}

}  // namespace Experimental
}  // namespace KokkosBlas

#endif
