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

#ifndef KOKKOSBLAS1_TEAM_ABS_SPEC_HPP_
#define KOKKOSBLAS1_TEAM_ABS_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Experimental {
namespace Impl {

template <class RV, class XV>
struct team_abs_tpl_spec_avail {
  constexpr static bool value = false;
};

// Unification and Specialization layer
template <class TeamType, class RV, class XV, bool tpl_spec_avail = team_abs_tpl_spec_avail<RV, XV>::value>
struct TeamAbs {
  typedef Kokkos::ArithTraits<typename XV::non_const_value_type> ATS;

  static KOKKOS_INLINE_FUNCTION void team_abs(const TeamType& team, const RV& R, const XV& X);
};

template <class TeamType, class RV, class XV>
struct TeamAbs<TeamType, RV, XV, false> {
  typedef Kokkos::ArithTraits<typename XV::non_const_value_type> ATS;

  static KOKKOS_INLINE_FUNCTION void team_abs(const TeamType& team, const RV& R, const XV& X) {
    int N = X.extent(0);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) { R(i) = ATS::abs(X(i)); });
  }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosBlas

#endif
