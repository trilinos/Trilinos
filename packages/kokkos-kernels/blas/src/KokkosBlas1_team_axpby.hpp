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
      team, a, x, Kokkos::ArithTraits<typename YVector::non_const_value_type>::one(), y);
}

}  // namespace Experimental
}  // namespace KokkosBlas

#endif
