/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSBLAS1_TEAM_NRM2_SPEC_HPP_
#define KOKKOSBLAS1_TEAM_NRM2_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Experimental {
namespace Impl {


template<class XV>
struct team_nrm2_tpl_spec_avail {
  constexpr static bool value = false;
};

// Unification and Specialization layer
template<class TeamType, class XV, bool tpl_spec_avail = team_nrm2_tpl_spec_avail<XV>::value>
struct TeamNrm2 {
  typedef typename Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type>::mag_type mag_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;
  
  static KOKKOS_INLINE_FUNCTION mag_type team_nrm2 (const TeamType& team, const XV& X);
};

template<class TeamType, class XV>
struct TeamNrm2<TeamType, XV, false> {
  typedef typename Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type>::mag_type mag_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<typename XV::non_const_value_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;

  static KOKKOS_INLINE_FUNCTION mag_type team_nrm2 (const TeamType& team, const XV& X) {
    mag_type result;
    int N = X.extent(0);
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,N), [&] (const int& i, mag_type& val) {
      const typename IPT::mag_type tmp = IPT::norm (X(i));
      val += tmp * tmp;
    },result);
    result = AT::sqrt(result);
    return result;
  }
};

}
}
}

#endif
