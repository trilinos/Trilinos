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

#ifndef KOKKOSBLAS1_TEAM_SCAL_SPEC_HPP_
#define KOKKOSBLAS1_TEAM_SCAL_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Experimental {
namespace Impl {

template <class RV, class XV>
struct team_scal_tpl_spec_avail {
  constexpr static bool value = false;
};

// Unification and Specialization layer
template <class TeamType, class RV, class XV, bool tpl_spec_avail = team_scal_tpl_spec_avail<RV, XV>::value>
struct TeamScal {
  static KOKKOS_INLINE_FUNCTION void team_scal(const TeamType& team, const RV& R,
                                               const typename XV::non_const_value_type& a, const XV& X);
};

template <class TeamType, class RV, class XV>
struct TeamScal<TeamType, RV, XV, false> {
  static KOKKOS_INLINE_FUNCTION void team_scal(const TeamType& team, const RV& R,
                                               const typename XV::non_const_value_type& a, const XV& X) {
    const int N = X.extent(0);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) { R(i) = a * X(i); });
  }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosBlas

#endif
