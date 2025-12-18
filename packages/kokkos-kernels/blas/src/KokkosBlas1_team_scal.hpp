// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_TEAM_SCAL_HPP_
#define KOKKOSBLAS1_TEAM_SCAL_HPP_

#include <KokkosBlas1_team_scal_spec.hpp>

namespace KokkosBlas {
namespace Experimental {

template <class TeamType, class RVector, class XVector>
void KOKKOS_INLINE_FUNCTION scal(const TeamType& team, const RVector& r,
                                 const typename XVector::non_const_value_type& a, const XVector& x) {
  return Impl::TeamScal<TeamType, RVector, XVector>::team_scal(team, r, a, x);
}

}  // namespace Experimental
}  // namespace KokkosBlas

#endif
