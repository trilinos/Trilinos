// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_TEAM_AXPBY_HPP_
#define KOKKOSBLAS1_TEAM_AXPBY_HPP_

#include <KokkosBlas1_team_axpby_spec.hpp>

namespace KokkosBlas {
namespace Experimental {

template <class TeamType, class XVector, class YVector>
void KOKKOS_INLINE_FUNCTION axpby(const TeamType& team, const typename XVector::non_const_value_type& a,
                                  const XVector& x, const typename YVector::non_const_value_type& b, const YVector& y) {
  return Impl::TeamAXPBY<TeamType, XVector, YVector>::team_axpby(team, a, x, b, y);
}

template <class TeamType, class XVector, class YVector>
void KOKKOS_INLINE_FUNCTION axpy(const TeamType& team, const typename XVector::non_const_value_type& a,
                                 const XVector& x, const YVector& y) {
  KokkosBlas::Experimental::axpby<TeamType, XVector, YVector>(
      team, a, x, KokkosKernels::ArithTraits<typename YVector::non_const_value_type>::one(), y);
}

}  // namespace Experimental
}  // namespace KokkosBlas

#endif
