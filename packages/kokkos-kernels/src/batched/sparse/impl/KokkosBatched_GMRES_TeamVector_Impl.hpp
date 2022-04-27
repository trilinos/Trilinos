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
#ifndef __KOKKOSBATCHED_GMRES_TEAMVECTOR_IMPL_HPP__
#define __KOKKOSBATCHED_GMRES_TEAMVECTOR_IMPL_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Axpy.hpp"
#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBatched_Dot.hpp"
#include "KokkosBatched_Spmv.hpp"
#include "KokkosBatched_Xpay.hpp"
#include "KokkosBatched_Givens_Serial_Internal.hpp"
#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_Identity.hpp"

namespace KokkosBatched {

///
/// TeamVector GMRES
///   Two nested parallel_for with both TeamVectorRange and ThreadVectorRange
///   (or one with TeamVectorRange) are used inside.
///

template <typename MemberType>
struct TeamVectorGMRES {
  template <typename OperatorType, typename VectorViewType,
            typename PrecOperatorType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType& member, const OperatorType& A, const VectorViewType& _B,
      const VectorViewType& _X, const PrecOperatorType& P,
      const KrylovHandle<typename VectorViewType::non_const_value_type>&
          handle) {
    typedef int OrdinalType;
    typedef typename Kokkos::Details::ArithTraits<
        typename VectorViewType::non_const_value_type>::mag_type MagnitudeType;
    typedef Kokkos::Details::ArithTraits<MagnitudeType> ATM;

    using ScratchPadNormViewType = Kokkos::View<
        MagnitudeType*,
        typename VectorViewType::execution_space::scratch_memory_space>;
    using ScratchPadVectorViewType = Kokkos::View<
        typename VectorViewType::non_const_value_type**,
        typename VectorViewType::array_layout,
        typename VectorViewType::execution_space::scratch_memory_space>;
    using ScratchPadMultiVectorViewType = Kokkos::View<
        typename VectorViewType::non_const_value_type***,
        typename VectorViewType::array_layout,
        typename VectorViewType::execution_space::scratch_memory_space>;
    using TeamVectorCopy1D = TeamVectorCopy<MemberType, Trans::NoTranspose, 1>;

    const OrdinalType numMatrices = _X.extent(0);
    const OrdinalType numRows     = _X.extent(1);

    size_t maximum_iteration = handle.get_max_iteration() < numRows
                                   ? handle.get_max_iteration()
                                   : numRows;
    const MagnitudeType tolerance     = handle.get_tolerance();
    const MagnitudeType max_tolerance = 0.;

    ScratchPadMultiVectorViewType V(member.team_scratch(1), numMatrices,
                                    maximum_iteration + 1, numRows);
    ScratchPadMultiVectorViewType H(member.team_scratch(1), numMatrices,
                                    maximum_iteration + 1, maximum_iteration);
    ScratchPadMultiVectorViewType Givens(member.team_scratch(1), numMatrices,
                                         maximum_iteration, 2);
    ScratchPadVectorViewType G(member.team_scratch(1), numMatrices,
                               maximum_iteration + 1);

    ScratchPadVectorViewType W(member.team_scratch(0), numMatrices, numRows);
    ScratchPadVectorViewType Q(member.team_scratch(0), numMatrices, numRows);
    ScratchPadVectorViewType R(member.team_scratch(0), numMatrices, numRows);
    ScratchPadVectorViewType X(member.team_scratch(0), numMatrices, numRows);

    ScratchPadNormViewType beta(member.team_scratch(0), numMatrices);
    ScratchPadNormViewType mask(member.team_scratch(0), numMatrices);
    ScratchPadNormViewType tmp(member.team_scratch(0), numMatrices);

    TeamVectorCopy<MemberType>::invoke(member, _X, X);
    // Deep copy of b into r_0:
    TeamVectorCopy<MemberType>::invoke(member, _B, R);

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, numMatrices),
                         [&](const OrdinalType& i) { mask(i) = 1.; });

    // r_0 := b - A x_0
    member.team_barrier();
    A.template apply<MemberType, ScratchPadVectorViewType,
                     ScratchPadVectorViewType, Trans::NoTranspose,
                     Mode::TeamVector>(member, X, R, -1, 1);
    member.team_barrier();

    P.template apply<MemberType, ScratchPadVectorViewType,
                     ScratchPadVectorViewType, Trans::NoTranspose,
                     Mode::TeamVector, 1>(member, R, R);
    member.team_barrier();

    TeamVectorDot<MemberType>::invoke(member, R, R, beta);
    member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, numMatrices),
                         [&](const OrdinalType& i) {
                           beta(i) = ATM::sqrt(beta(i));
                           G(i, 0) = beta(i) > max_tolerance ? beta(i) : 0.;
                           tmp(i) = beta(i) > max_tolerance ? 1. / beta(i) : 0.;
                         });

    member.team_barrier();  // Finish writing to tmp

    Kokkos::parallel_for(
        Kokkos::TeamVectorRange(member, 0, numMatrices * numRows),
        [&](const OrdinalType& iTemp) {
          OrdinalType iRow, iMatrix;
          getIndices<OrdinalType, typename VectorViewType::array_layout>(
              iTemp, numRows, numMatrices, iRow, iMatrix);
          V(iMatrix, 0, iRow) = R(iMatrix, iRow) * tmp(iMatrix);
        });

    int status = 1;
    // int number_not_converged = 0;

    for (size_t j = 0; j < maximum_iteration; ++j) {
      member.team_barrier();  // Finish writing to V
      // q := A p_j
      auto V_j = Kokkos::subview(V, Kokkos::ALL, j, Kokkos::ALL);

      A.template apply<MemberType, ScratchPadVectorViewType,
                       ScratchPadVectorViewType, Trans::NoTranspose,
                       Mode::TeamVector>(member, V_j, W);
      member.team_barrier();
      P.template apply<MemberType, ScratchPadVectorViewType,
                       ScratchPadVectorViewType, Trans::NoTranspose,
                       Mode::TeamVector, 1>(member, W, W);

      for (size_t i = 0; i < j + 1; ++i) {
        member.team_barrier();  // Finish writing to W
        auto V_i = Kokkos::subview(V, Kokkos::ALL, i, Kokkos::ALL);
        TeamVectorDot<MemberType>::invoke(member, W, V_i, tmp);
        member.team_barrier();
        TeamVectorCopy1D::invoke(member, tmp,
                                 Kokkos::subview(H, Kokkos::ALL, i, j));

        member.team_barrier();  // Don't start modifying tmp until copy above
                                // finishes
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(member, 0, numMatrices),
            [&](const OrdinalType& ii) { tmp(ii) = -tmp(ii); });

        member.team_barrier();  // Finish writing to tmp

        TeamVectorAxpy<MemberType>::invoke(member, tmp, V_i, W);
      }

      member.team_barrier();  // Finish writing to W
      TeamVectorDot<MemberType>::invoke(member, W, W, tmp);
      member.team_barrier();
      Kokkos::parallel_for(
          Kokkos::TeamVectorRange(member, 0, numMatrices),
          [&](const OrdinalType& i) {
            H(i, j + 1, j) = ATM::sqrt(tmp(i));
            tmp(i) = H(i, j + 1, j) > max_tolerance ? 1. / H(i, j + 1, j) : 0.;
          });
      member.team_barrier();
      Kokkos::parallel_for(
          Kokkos::TeamVectorRange(member, 0, numMatrices * numRows),
          [&](const OrdinalType& iTemp) {
            OrdinalType iRow, iMatrix;
            getIndices<OrdinalType, typename VectorViewType::array_layout>(
                iTemp, numRows, numMatrices, iRow, iMatrix);
            V(iMatrix, j + 1, iRow) = W(iMatrix, iRow) * tmp(iMatrix);
          });

      Kokkos::parallel_for(
          Kokkos::TeamVectorRange(member, 0, numMatrices),
          [&](const OrdinalType& l) {
            // Apply the previous Givens rotations:
            auto H_j = Kokkos::subview(H, l, Kokkos::ALL, j);

            if (mask(l) == 1.) {
              for (size_t i = 0; i < j; ++i) {
                auto tmp1 =
                    Givens(l, i, 0) * H_j(i) + Givens(l, i, 1) * H_j(i + 1);
                auto tmp2 =
                    -Givens(l, i, 1) * H_j(i) + Givens(l, i, 0) * H_j(i + 1);
                H_j(i)     = tmp1;
                H_j(i + 1) = tmp2;
              }

              // Compute the new Givens rotation:
              Kokkos::pair<typename VectorViewType::non_const_value_type,
                           typename VectorViewType::non_const_value_type>
                  G_new(1, 0);
              typename VectorViewType::non_const_value_type alpha = 0;
              SerialGivensInternal::invoke(H_j(j), H_j(j + 1), &G_new, &alpha);

              Givens(l, j, 0) = G_new.first;
              Givens(l, j, 1) = G_new.second;

              // Apply the new Givens rotation:
              auto tmp1 =
                  Givens(l, j, 0) * H_j(j) + Givens(l, j, 1) * H_j(j + 1);
              auto tmp2 =
                  -Givens(l, j, 1) * H_j(j) + Givens(l, j, 0) * H_j(j + 1);
              H_j(j)     = tmp1;
              H_j(j + 1) = tmp2;

              G(l, j + 1) = -Givens(l, j, 1) * G(l, j);
              G(l, j) *= Givens(l, j, 0);
            } else {
              H_j(j)      = 1.;
              G(l, j + 1) = 0.;
            }

            if (mask(l) == 1. &&
                Kokkos::ArithTraits<double>::abs(G(l, j + 1)) / beta(l) <
                    tolerance) {
              mask(l)     = 0.;
              G(l, j + 1) = 0.;
            }
          });
    }

    member.team_barrier();  // Finish writing to G

    Kokkos::parallel_for(
        Kokkos::TeamVectorRange(member, 0, numMatrices),
        [&](const OrdinalType& l) {
          SerialTrsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit,
                     Algo::Trsm::Unblocked>::template invoke(1,
                                                             Kokkos::subview(
                                                                 H, l,
                                                                 Kokkos::ALL,
                                                                 Kokkos::ALL),
                                                             Kokkos::subview(
                                                                 G, l,
                                                                 Kokkos::ALL));
        });

    member.team_barrier();  // Finish writing to G

    for (size_t j = 0; j < maximum_iteration; ++j) {
      TeamVectorAxpy<MemberType>::invoke(
          member, Kokkos::subview(G, Kokkos::ALL, j),
          Kokkos::subview(V, Kokkos::ALL, j, Kokkos::ALL), X);
      member.team_barrier();  // Finish writing to X
    }

    TeamVectorCopy<MemberType>::invoke(member, X, _X);
    return status;
  }

  template <typename OperatorType, typename VectorViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType& member, const OperatorType& A, const VectorViewType& _B,
      const VectorViewType& _X,
      const KrylovHandle<typename VectorViewType::non_const_value_type>&
          handle) {
    Identity P;
    return invoke<OperatorType, VectorViewType, Identity>(member, A, _B, _X, P,
                                                          handle);
  }
};
}  // namespace KokkosBatched

#endif
