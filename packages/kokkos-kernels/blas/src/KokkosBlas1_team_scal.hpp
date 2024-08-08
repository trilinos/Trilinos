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

#ifndef KOKKOSBLAS1_TEAM_SCAL_HPP_
#define KOKKOSBLAS1_TEAM_SCAL_HPP_

#include <KokkosBlas1_team_scal_spec.hpp>

namespace KokkosBlas {
namespace Experimental {

template <class TeamType, class RVector, class XVector>
void KOKKOS_INLINE_FUNCTION scal(const TeamType& team, const RVector& r,
                                 const typename XVector::non_const_value_type& a, const XVector& x) {
  return Impl::TeamScal<TeamType, RVector, XVector>::team_scal(team, r, a, x);
}

}  // namespace Experimental
}  // namespace KokkosBlas

#endif
