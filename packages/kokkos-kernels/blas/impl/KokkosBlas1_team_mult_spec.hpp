// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_TEAM_MULT_SPEC_HPP_
#define KOKKOSBLAS1_TEAM_MULT_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Experimental {
namespace Impl {

template <class YV, class AV, class XV>
struct team_mult_tpl_spec_avail {
  constexpr static bool value = false;
};

// Unification and Specialization layer
template <class TeamType, class YVector, class AVector, class XVector,
          bool tpl_spec_avail = team_mult_tpl_spec_avail<YVector, AVector, XVector>::value>
struct TeamMult {
  static KOKKOS_INLINE_FUNCTION void team_mult(const TeamType& team,
                                               const typename YVector::non_const_value_type& gamma, const YVector& y,
                                               const typename AVector::non_const_value_type& alpha, const AVector& a,
                                               const XVector& x);
};

template <class TeamType, class YVector, class AVector, class XVector>
struct TeamMult<TeamType, YVector, AVector, XVector, false> {
  static KOKKOS_INLINE_FUNCTION void team_mult(const TeamType& team,
                                               const typename YVector::non_const_value_type& gamma, const YVector& y,
                                               const typename AVector::non_const_value_type& alpha, const AVector& a,
                                               const XVector& x) {
    const int N = x.extent(0);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N),
                         [&](const int& i) { y(i) = gamma * y(i) + alpha * a(i) * x(i); });
  }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosBlas

#endif
