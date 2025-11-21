// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_TEAM_NRM2_HPP_
#define KOKKOSBLAS1_TEAM_NRM2_HPP_

#include <KokkosBlas1_team_nrm2_spec.hpp>

namespace KokkosBlas {
namespace Experimental {

template <class TeamType, class XVector>
typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::mag_type
    KOKKOS_INLINE_FUNCTION
    nrm2(const TeamType& team, const XVector& x) {
  return Impl::TeamNrm2<TeamType, XVector>::team_nrm2(team, x);
}

}  // namespace Experimental
}  // namespace KokkosBlas

#endif
