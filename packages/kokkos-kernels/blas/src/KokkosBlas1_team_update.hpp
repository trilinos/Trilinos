// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_TEAM_UPDATE_HPP_
#define KOKKOSBLAS1_TEAM_UPDATE_HPP_

#include <KokkosBlas1_team_update_spec.hpp>

namespace KokkosBlas {
namespace Experimental {

template <class TeamType, class XVector, class YVector, class ZVector>
void KOKKOS_INLINE_FUNCTION update(const TeamType& team, const typename XVector::non_const_value_type& alpha,
                                   const XVector& x, const typename YVector::non_const_value_type& beta,
                                   const YVector& y, const typename ZVector::non_const_value_type& gamma,
                                   const ZVector& z) {
  return Impl::TeamUpdate<TeamType, XVector, YVector, ZVector>::team_update(team, alpha, x, beta, y, gamma, z);
}

}  // namespace Experimental
}  // namespace KokkosBlas

#endif
