// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_TEAM_MULT_HPP_
#define KOKKOSBLAS1_TEAM_MULT_HPP_

#include <KokkosBlas1_team_mult_spec.hpp>

namespace KokkosBlas {
namespace Experimental {

template <class TeamType, class YVector, class AVector, class XVector>
void KOKKOS_INLINE_FUNCTION mult(const TeamType& team, const typename YVector::non_const_value_type& gamma,
                                 const YVector& y, const typename AVector::non_const_value_type& alpha,
                                 const AVector& a, const XVector& x) {
  return Impl::TeamMult<TeamType, YVector, AVector, XVector>::team_mult(team, gamma, y, alpha, a, x);
}

}  // namespace Experimental
}  // namespace KokkosBlas

#endif
