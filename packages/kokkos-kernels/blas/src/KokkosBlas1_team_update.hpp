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
