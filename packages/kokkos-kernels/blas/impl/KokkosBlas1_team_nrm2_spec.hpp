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

#ifndef KOKKOSBLAS1_TEAM_NRM2_SPEC_HPP_
#define KOKKOSBLAS1_TEAM_NRM2_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Experimental {
namespace Impl {

template <class XV>
struct team_nrm2_tpl_spec_avail {
  constexpr static bool value = false;
};

// Unification and Specialization layer
template <class TeamType, class XV, bool tpl_spec_avail = team_nrm2_tpl_spec_avail<XV>::value>
struct TeamNrm2 {
  typedef typename Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type>::mag_type mag_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type> IPT;
  typedef Kokkos::ArithTraits<typename IPT::mag_type> AT;

  static KOKKOS_INLINE_FUNCTION mag_type team_nrm2(const TeamType& team, const XV& X);
};

template <class TeamType, class XV>
struct TeamNrm2<TeamType, XV, false> {
  typedef typename Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type>::mag_type mag_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type> IPT;
  typedef Kokkos::ArithTraits<typename IPT::mag_type> AT;

  static KOKKOS_INLINE_FUNCTION mag_type team_nrm2(const TeamType& team, const XV& X) {
    mag_type result = 0.0;  // Kokkos::ArithTraits<mag_type>zero();
    int N           = X.extent(0);
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, N),
        [&](const int& i, mag_type& val) {
          const typename IPT::mag_type tmp = IPT::norm(X(i));
          val += tmp * tmp;
        },
        result);
    result = AT::sqrt(result);
    return result;
  }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace KokkosBlas

#endif
