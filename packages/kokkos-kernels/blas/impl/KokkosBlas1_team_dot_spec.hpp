// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_TEAM_DOT_SPEC_HPP_
#define KOKKOSBLAS1_TEAM_DOT_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Experimental {
namespace Impl {

template <class XV, class YV>
struct team_dot_tpl_spec_avail {
  constexpr static bool value = false;
};

// Unification and Specialization layer
template <class TeamType, class XV, class YV, bool tpl_spec_avail = team_dot_tpl_spec_avail<XV, YV>::value>
struct TeamDot {
  typedef Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type> IPT;
  typedef typename IPT::dot_type dot_type;

  static KOKKOS_INLINE_FUNCTION dot_type team_dot(const TeamType& team, const XV& X, const YV& Y);
};

template <class TeamType, class XV, class YV>
struct TeamDot<TeamType, XV, YV, false> {
  typedef Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type> IPT;
  typedef typename IPT::dot_type dot_type;

  static KOKKOS_INLINE_FUNCTION dot_type team_dot(const TeamType& team, const XV& X, const YV& Y) {
    dot_type result = 0.0;  // KokkosKernels::ArithTraits<dot_type>zero();
    int N           = X.extent(0);
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, N),
        [&](const int& i, dot_type& val) {
          val += IPT::dot(X(i), Y(i));  // X(i) * Y(i)
        },
        result);
    return result;
  }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosBlas

#endif
