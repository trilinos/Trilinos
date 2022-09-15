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
#ifndef __KOKKOSBATCHED_GMRES_TEAM_IMPL_HPP__
#define __KOKKOSBATCHED_GMRES_TEAM_IMPL_HPP__

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
#include "KokkosBatched_Gemv_Decl.hpp"

namespace KokkosBatched {

///
/// Team GMRES
///   A nested parallel_for with TeamThreadRange is used.
///

template <typename MemberType>
template <typename OperatorType, typename VectorViewType,
          typename PrecOperatorType, typename KrylovHandleType>
KOKKOS_INLINE_FUNCTION int TeamGMRES<MemberType>::invoke(
    const MemberType& member, const OperatorType& A, const VectorViewType& _B,
    const VectorViewType& _X, const PrecOperatorType& P,
    const KrylovHandleType& handle) {
  typedef int OrdinalType;
  typedef typename Kokkos::Details::ArithTraits<
      typename VectorViewType::non_const_value_type>::mag_type MagnitudeType;
  typedef Kokkos::Details::ArithTraits<MagnitudeType> ATM;

  using ScratchPadVectorViewType = Kokkos::View<
      typename VectorViewType::non_const_value_type**,
      typename VectorViewType::array_layout,
      typename VectorViewType::execution_space::scratch_memory_space>;
  using TeamCopy1D = TeamCopy<MemberType, Trans::NoTranspose, 1>;

  const OrdinalType numMatrices = _X.extent(0);
  const OrdinalType numRows     = _X.extent(1);

  size_t maximum_iteration = handle.get_max_iteration() < numRows
                                 ? handle.get_max_iteration()
                                 : numRows;
  const MagnitudeType tolerance     = handle.get_tolerance();
  const MagnitudeType max_tolerance = handle.get_max_tolerance();

  int n_V      = numRows;
  int n_H      = maximum_iteration + 1;
  int n_Givens = 2;

  int offset_V      = 0;
  int offset_H      = offset_V + n_V;
  int offset_Givens = offset_H + n_H;

  const int first_matrix = handle.first_index(member.league_rank());
  const int last_matrix  = handle.last_index(member.league_rank());

  auto V_view = Kokkos::subview(
      handle.Arnoldi_view, Kokkos::make_pair(first_matrix, last_matrix),
      Kokkos::ALL, Kokkos::make_pair(offset_V, offset_V + n_V));
  auto H_view = Kokkos::subview(
      handle.Arnoldi_view, Kokkos::make_pair(first_matrix, last_matrix),
      Kokkos::ALL, Kokkos::make_pair(offset_H, offset_H + n_H));
  auto Givens_view = Kokkos::subview(
      handle.Arnoldi_view, Kokkos::make_pair(first_matrix, last_matrix),
      Kokkos::ALL, Kokkos::make_pair(offset_Givens, offset_Givens + n_Givens));

  int n_G    = maximum_iteration + 1;
  int n_W    = numRows;
  int n_X    = numRows;
  int n_mask = 1;
  int n_tmp  = 1;

  int offset_G    = 0;
  int offset_W    = offset_G + n_G;
  int offset_X    = offset_W + n_W;
  int offset_mask = offset_X + n_X;
  int offset_tmp  = offset_mask + n_mask;

  ScratchPadVectorViewType tmp_2D(
      member.team_scratch(handle.get_scratch_pad_level()), numMatrices,
      n_G + n_W + n_X + n_mask + n_tmp);

  auto G    = Kokkos::subview(tmp_2D, Kokkos::ALL,
                           Kokkos::make_pair(offset_G, offset_G + n_G));
  auto W    = Kokkos::subview(tmp_2D, Kokkos::ALL,
                           Kokkos::make_pair(offset_W, offset_W + n_W));
  auto X    = Kokkos::subview(tmp_2D, Kokkos::ALL,
                           Kokkos::make_pair(offset_X, offset_X + n_X));
  auto mask = Kokkos::subview(tmp_2D, Kokkos::ALL, offset_mask);
  auto tmp  = Kokkos::subview(tmp_2D, Kokkos::ALL, offset_tmp);

  TeamCopy<MemberType>::invoke(member, _X, X);
  // Deep copy of b into r_0:
  TeamCopy<MemberType>::invoke(member, _B, W);

  // r_0 := b - A x_0
  member.team_barrier();
  A.template apply<Trans::NoTranspose, Mode::Team>(member, X, W, -1, 1);
  member.team_barrier();

  P.template apply<Trans::NoTranspose, Mode::Team, 1>(member, W, W);
  member.team_barrier();

  TeamDot<MemberType>::invoke(member, W, W, tmp);
  member.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, numMatrices),
                       [&](const OrdinalType& i) {
                         tmp(i) = ATM::sqrt(tmp(i));
                         handle.set_norm(member.league_rank(), i, 0, tmp(i));
                         if (tmp(i) > max_tolerance) {
                           mask(i) = 1;
                           G(i, 0) = tmp(i);
                           tmp(i)  = 1. / tmp(i);
                         } else {
                           handle.set_iteration(member.league_rank(), i, 0);
                           mask(i) = 0;
                           G(i, 0) = 0.;
                           tmp(i)  = 0.;
                         }
                       });

  member.team_barrier();  // Finish writing to tmp

  auto V_0 = Kokkos::subview(V_view, Kokkos::ALL, 0, Kokkos::ALL);
  Kokkos::parallel_for(
      Kokkos::TeamThreadRange(member, 0, numMatrices * numRows),
      [&](const OrdinalType& iTemp) {
        OrdinalType iRow, iMatrix;
        getIndices<OrdinalType, typename VectorViewType::array_layout>(
            iTemp, numRows, numMatrices, iRow, iMatrix);
        V_0(iMatrix, iRow) = W(iMatrix, iRow) * tmp(iMatrix);
      });
  int status = 1;
  // int number_not_converged = 0;

  for (size_t j = 0; j < maximum_iteration; ++j) {
    member.team_barrier();  // Finish writing to V
    // q := A p_j
    auto V_j = Kokkos::subview(V_view, Kokkos::ALL, j, Kokkos::ALL);

    A.template apply<Trans::NoTranspose, Mode::Team>(member, V_j, W);
    member.team_barrier();

    P.template apply<Trans::NoTranspose, Mode::Team, 1>(member, W, W);
    member.team_barrier();

    if (handle.get_ortho_strategy() == 0) {
      auto V_old = Kokkos::subview(
          V_view, Kokkos::ALL, Kokkos::make_pair(0, (int)j + 1), Kokkos::ALL);
      auto H_old = Kokkos::subview(H_view, Kokkos::ALL, j,
                                   Kokkos::make_pair(0, (int)j + 1));
      // Inner products
      TeamGemv<MemberType, Trans::NoTranspose, Algo::Gemv::Unblocked>::invoke(
          member, 1, V_old, W, 0, H_old);
      member.team_barrier();

      // Update
      TeamGemv<MemberType, Trans::Transpose, Algo::Gemv::Unblocked>::invoke(
          member, -1, V_old, H_old, 1, W);
      member.team_barrier();  // Finish writing to W
    }
    if (handle.get_ortho_strategy() == 1) {
      for (size_t i = 0; i < j + 1; ++i) {
        auto V_i = Kokkos::subview(V_view, Kokkos::ALL, i, Kokkos::ALL);
        TeamDot<MemberType>::invoke(member, W, V_i, tmp);
        member.team_barrier();
        TeamCopy1D::invoke(member, tmp,
                           Kokkos::subview(H_view, Kokkos::ALL, j, i));
        member.team_barrier();
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(member, 0, numMatrices),
            [&](const OrdinalType& ii) { tmp(ii) = -tmp(ii); });

        member.team_barrier();  // Finish writing to tmp

        TeamAxpy<MemberType>::invoke(member, tmp, V_i, W);
        member.team_barrier();  // Finish writing to W
      }
    }

    TeamDot<MemberType>::invoke(member, W, W, tmp);
    member.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, numMatrices),
                         [&](const OrdinalType& i) {
                           H_view(i, j, j + 1) = ATM::sqrt(tmp(i));
                           tmp(i) = H_view(i, j, j + 1) > max_tolerance
                                        ? 1. / H_view(i, j, j + 1)
                                        : 0.;
                         });
    member.team_barrier();
    if (j + 1 < maximum_iteration) {
      auto V_n = Kokkos::subview(V_view, Kokkos::ALL, j + 1, Kokkos::ALL);
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, 0, numMatrices * numRows),
          [&](const OrdinalType& iTemp) {
            OrdinalType iRow, iMatrix;
            getIndices<OrdinalType, typename VectorViewType::array_layout>(
                iTemp, numRows, numMatrices, iRow, iMatrix);
            V_n(iMatrix, iRow) = W(iMatrix, iRow) * tmp(iMatrix);
          });
      member.team_barrier();
    }

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(member, 0, numMatrices),
        [&](const OrdinalType& l) {
          // Apply the previous Givens rotations:
          auto H_j        = Kokkos::subview(H_view, l, j, Kokkos::ALL);
          auto Givens_0_l = Kokkos::subview(Givens_view, l, Kokkos::ALL, 0);
          auto Givens_1_l = Kokkos::subview(Givens_view, l, Kokkos::ALL, 1);

          if (mask(l) == 1.) {
            for (size_t i = 0; i < j; ++i) {
              auto tmp1  = Givens_0_l(i) * H_j(i) + Givens_1_l(i) * H_j(i + 1);
              auto tmp2  = -Givens_1_l(i) * H_j(i) + Givens_0_l(i) * H_j(i + 1);
              H_j(i)     = tmp1;
              H_j(i + 1) = tmp2;
            }

            // Compute the new Givens rotation:
            Kokkos::pair<typename VectorViewType::non_const_value_type,
                         typename VectorViewType::non_const_value_type>
                G_new(1, 0);
            typename VectorViewType::non_const_value_type alpha = 0;
            SerialGivensInternal::invoke(H_j(j), H_j(j + 1), &G_new, &alpha);

            Givens_0_l(j) = G_new.first;
            Givens_1_l(j) = G_new.second;

            // Apply the new Givens rotation:
            auto tmp1  = Givens_0_l(j) * H_j(j) + Givens_1_l(j) * H_j(j + 1);
            auto tmp2  = -Givens_1_l(j) * H_j(j) + Givens_0_l(j) * H_j(j + 1);
            H_j(j)     = tmp1;
            H_j(j + 1) = tmp2;

            G(l, j + 1) = -Givens_1_l(j) * G(l, j);
            G(l, j) *= Givens_0_l(j);
          } else {
            H_j(j)      = 1.;
            G(l, j + 1) = 0.;
          }

          auto res_norm =
              Kokkos::ArithTraits<double>::abs(G(l, j + 1)) / G(l, 0);

          handle.set_norm(member.league_rank(), l, j + 1, res_norm);

          if (mask(l) == 1. && res_norm < tolerance) {
            mask(l)     = 0.;
            G(l, j + 1) = 0.;
            handle.set_iteration(member.league_rank(), l, j + 1);
          }
        });
    member.team_barrier();

    bool all_converged = true;
    for (OrdinalType l = 0; l < numMatrices; ++l)
      all_converged = (all_converged && mask(l) == 0.);
    if (all_converged) {
      maximum_iteration = j + 1;
      break;
    }
  }

  member.team_barrier();  // Finish writing to G

  auto first_indices = Kokkos::make_pair(0, (int)maximum_iteration);

  Kokkos::parallel_for(
      Kokkos::TeamVectorRange(member, 0, numMatrices),
      [&](const OrdinalType& l) {
        auto A_l = Kokkos::subview(H_view, l, first_indices, first_indices);
        auto B_l = Kokkos::subview(G, l, first_indices);

        SerialTrsm<Side::Left, Uplo::Lower, Trans::Transpose, Diag::NonUnit,
                   Algo::Trsm::Unblocked>::invoke(1, A_l, B_l);
      });

  member.team_barrier();  // Finish writing to G

  if (handle.get_ortho_strategy() == 0) {
    TeamGemv<MemberType, Trans::Transpose, Algo::Gemv::Unblocked>::invoke(
        member, 1,
        Kokkos::subview(V_view, Kokkos::ALL, first_indices, Kokkos::ALL),
        Kokkos::subview(G, Kokkos::ALL, first_indices), 1, X);
    member.team_barrier();  // Finish writing to X
  }
  if (handle.get_ortho_strategy() == 1) {
    for (size_t j = 0; j < maximum_iteration; ++j) {
      TeamAxpy<MemberType>::invoke(
          member, Kokkos::subview(G, Kokkos::ALL, j),
          Kokkos::subview(V_view, Kokkos::ALL, j, Kokkos::ALL), X);
      member.team_barrier();  // Finish writing to X
    }
  }

  TeamCopy<MemberType>::invoke(member, X, _X);

  member.team_barrier();

  if (handle.get_compute_last_residual()) {
    TeamCopy<MemberType>::invoke(member, _B, W);
    member.team_barrier();
    A.template apply<Trans::NoTranspose, Mode::Team>(member, X, W, -1, 1);
    member.team_barrier();
    P.template apply<Trans::NoTranspose, Mode::Team, 1>(member, W, W);
    member.team_barrier();
    TeamDot<MemberType>::invoke(member, W, W, tmp);
    member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, numMatrices),
                         [&](const OrdinalType& i) {
                           tmp(i) = ATM::sqrt(tmp(i));
                           handle.set_last_norm(member.league_rank(), i,
                                                tmp(i));
                         });
  }
  return status;
}

template <typename MemberType>
template <typename OperatorType, typename VectorViewType,
          typename KrylovHandleType>
KOKKOS_INLINE_FUNCTION int TeamGMRES<MemberType>::invoke(
    const MemberType& member, const OperatorType& A, const VectorViewType& _B,
    const VectorViewType& _X, const KrylovHandleType& handle) {
  Identity P;
  return invoke<OperatorType, VectorViewType, Identity>(member, A, _B, _X, P,
                                                        handle);
}

}  // namespace KokkosBatched

#endif
