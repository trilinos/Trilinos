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
