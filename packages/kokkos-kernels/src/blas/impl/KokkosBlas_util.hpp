/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Brian Kelley (bmkelle@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_BLAS_UTIL_HPP
#define KOKKOS_BLAS_UTIL_HPP

namespace KokkosBlas {
namespace Impl {

// Helper to choose the work distribution for a TeamPolicy computing multiple
// reductions. Each team computes a partial reduction and atomically contributes
// to the final result.
//
// This was originally written for dot-based GEMM, but can also be applied to
// multivector dots/norms.

// Input params:
//  * length: size of each vector to reduce
//  * numReductions: number of reductions to compute
// Output params:
//  * teamsPerReduction: number of teams to use for each reduction
template <typename ExecSpace, typename size_type>
void multipleReductionWorkDistribution(size_type length,
                                       size_type numReductions,
                                       size_type& teamsPerDot) {
  constexpr size_type workPerTeam = 4096;  // Amount of work per team
  size_type appxNumTeams =
      (length * numReductions) / workPerTeam;  // Estimation for appxNumTeams

  // Adjust appxNumTeams in case it is too small or too large
  if (appxNumTeams < 1) appxNumTeams = 1;
  if (appxNumTeams > 1024) appxNumTeams = 1024;

  // If there are more reductions than the number of teams,
  // then set the number of teams to be number of reductions.
  // We don't want a team to contribute to more than one reduction.
  if (numReductions >= appxNumTeams) {
    teamsPerDot = 1;
  }
  // If there are more teams than reductions, each reduction can
  // potentially be performed by multiple teams. First, compute
  // teamsPerDot as an integer (take the floor, not ceiling), then,
  // compute actual number of teams by using this factor.
  else {
    teamsPerDot = appxNumTeams / numReductions;
  }
}

// Functor to apply sqrt() to each element of a 1D view.

template <class RV>
struct TakeSqrtFunctor {
  TakeSqrtFunctor(const RV& r_) : r(r_) {}

  KOKKOS_INLINE_FUNCTION void operator()(int i) const {
    r(i) = Kokkos::ArithTraits<typename RV::non_const_value_type>::sqrt(r(i));
  }

  RV r;
};

}  // namespace Impl
}  // namespace KokkosBlas

#endif
