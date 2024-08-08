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
