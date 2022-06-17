//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.4
//       Copyright (2021) National Technology & Engineering
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
#ifndef __KOKKOSBATCHED_CG_TEAM_IMPL_HPP__
#define __KOKKOSBATCHED_CG_TEAM_IMPL_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Axpy.hpp"
#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBatched_Dot.hpp"
#include "KokkosBatched_Spmv.hpp"
#include "KokkosBatched_Xpay.hpp"

namespace KokkosBatched {

///
/// Team CG
///   A nested parallel_for with TeamThreadRange is used.
///

template <typename MemberType>
struct TeamCG {
  template <typename OperatorType, typename VectorViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType& member, const OperatorType& A, const VectorViewType& _B,
      const VectorViewType& _X,
      const KrylovHandle<typename VectorViewType::non_const_value_type>&
          handle) {
    typedef int OrdinalType;
    typedef typename Kokkos::Details::ArithTraits<
        typename VectorViewType::non_const_value_type>::mag_type MagnitudeType;

    size_t maximum_iteration      = handle.get_max_iteration();
    const MagnitudeType tolerance = handle.get_tolerance();

    using ScratchPadNormViewType = Kokkos::View<
        MagnitudeType*,
        typename VectorViewType::execution_space::scratch_memory_space>;
    using ScratchPadVectorViewType = Kokkos::View<
        typename VectorViewType::non_const_value_type**,
        typename VectorViewType::array_layout,
        typename VectorViewType::execution_space::scratch_memory_space>;
    using TeamCopy1D = TeamCopy<MemberType, Trans::NoTranspose, 1>;

    const OrdinalType numMatrices = _X.extent(0);
    const OrdinalType numRows     = _X.extent(1);

    ScratchPadVectorViewType P(member.team_scratch(0), numMatrices, numRows);
    ScratchPadVectorViewType Q(member.team_scratch(0), numMatrices, numRows);
    ScratchPadVectorViewType R(member.team_scratch(0), numMatrices, numRows);
    ScratchPadVectorViewType X(member.team_scratch(0), numMatrices, numRows);

    ScratchPadNormViewType sqr_norm_0(member.team_scratch(0), numMatrices);
    ScratchPadNormViewType sqr_norm_j(member.team_scratch(0), numMatrices);
    ScratchPadNormViewType alpha(member.team_scratch(0), numMatrices);
    ScratchPadNormViewType mask(member.team_scratch(0), numMatrices);
    ScratchPadNormViewType tmp(member.team_scratch(0), numMatrices);

    TeamCopy<MemberType>::invoke(member, _X, X);
    // Deep copy of b into r_0:
    TeamCopy<MemberType>::invoke(member, _B, R);

    // r_0 := b - A x_0
    member.team_barrier();
    A.template apply<MemberType, ScratchPadVectorViewType,
                     ScratchPadVectorViewType, Trans::NoTranspose, Mode::Team>(
        member, X, R, -1, 1);
    member.team_barrier();

    // Deep copy of r_0 into p_0:
    TeamCopy<MemberType>::invoke(member, R, P);

    TeamDot<MemberType>::invoke(member, R, R, sqr_norm_0);
    member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, numMatrices),
                         [&](const OrdinalType& i) {
                           mask(i) =
                               sqr_norm_0(i) > tolerance * tolerance ? 1. : 0;
                         });

    TeamCopy1D::invoke(member, sqr_norm_0, sqr_norm_j);

    int status               = 1;
    int number_not_converged = 0;

    for (size_t j = 0; j < maximum_iteration; ++j) {
      // q := A p_j
      A.template apply<MemberType, ScratchPadVectorViewType,
                       ScratchPadVectorViewType, Trans::NoTranspose,
                       Mode::Team>(member, P, Q);
      member.team_barrier();

      TeamDot<MemberType>::invoke(member, P, Q, tmp);
      member.team_barrier();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, numMatrices),
                           [&](const OrdinalType& i) {
                             alpha(i) =
                                 mask(i) != 0. ? sqr_norm_j(i) / tmp(i) : 0.;
                           });
      member.team_barrier();

      // x_{j+1} := alpha p_j + x_j
      TeamAxpy<MemberType>::invoke(member, alpha, P, X);
      member.team_barrier();

      // r_{j+1} := - alpha q + r_j
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, numMatrices),
                           [&](const OrdinalType& i) { alpha(i) = -alpha(i); });
      member.team_barrier();

      TeamAxpy<MemberType>::invoke(member, alpha, Q, R);
      member.team_barrier();

      TeamDot<MemberType>::invoke(member, R, R, tmp);
      member.team_barrier();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, numMatrices),
                           [&](const OrdinalType& i) {
                             alpha(i) =
                                 mask(i) != 0. ? tmp(i) / sqr_norm_j(i) : 0.;
                           });

      TeamCopy1D::invoke(member, tmp, sqr_norm_j);

      // Relative convergence check:
      number_not_converged = 0;
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(member, 0, numMatrices),
          [&](const OrdinalType& i, int& lnumber_not_converged) {
            if (sqr_norm_j(i) / sqr_norm_0(i) > tolerance * tolerance)
              ++lnumber_not_converged;
            else
              mask(i) = 0.;
          },
          number_not_converged);

      member.team_barrier();

      if (number_not_converged == 0) {
        status = 0;
        break;
      }

      // p_{j+1} := alpha p_j + r_{j+1}
      TeamXpay<MemberType>::invoke(member, alpha, R, P);
      member.team_barrier();
    }

    TeamCopy<MemberType>::invoke(member, X, _X);
    return status;
  }
};
}  // namespace KokkosBatched

#endif
