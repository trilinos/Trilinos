// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_CG_TEAM_IMPL_HPP
#define KOKKOSBATCHED_CG_TEAM_IMPL_HPP

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
template <typename OperatorType, typename VectorViewType, typename KrylovHandle, typename TMPViewType,
          typename TMPNormViewType>
KOKKOS_INLINE_FUNCTION int TeamCG<MemberType>::invoke(const MemberType& member, const OperatorType& A,
                                                      const VectorViewType& B, const VectorViewType& X,
                                                      const KrylovHandle& handle, const TMPViewType& TMPView,
                                                      const TMPNormViewType& TMPNormView) {
  typedef int OrdinalType;
  typedef typename KokkosKernels::ArithTraits<typename VectorViewType::non_const_value_type>::mag_type MagnitudeType;

  size_t maximum_iteration      = handle.get_max_iteration();
  const MagnitudeType tolerance = handle.get_tolerance();

  using TeamCopy1D = TeamCopy<MemberType, Trans::NoTranspose, 1>;

  const OrdinalType numMatrices = X.extent(0);
  const OrdinalType numRows     = X.extent(1);

  int offset_P = 0;
  int offset_Q = offset_P + numRows;
  int offset_R = offset_Q + numRows;
  int offset_X = offset_R + numRows;

  auto P  = Kokkos::subview(TMPView, Kokkos::ALL, Kokkos::make_pair(offset_P, offset_P + numRows));
  auto Q  = Kokkos::subview(TMPView, Kokkos::ALL, Kokkos::make_pair(offset_Q, offset_Q + numRows));
  auto R  = Kokkos::subview(TMPView, Kokkos::ALL, Kokkos::make_pair(offset_R, offset_R + numRows));
  auto X_ = Kokkos::subview(TMPView, Kokkos::ALL, Kokkos::make_pair(offset_X, offset_X + numRows));

  auto sqr_norm_0 = Kokkos::subview(TMPNormView, Kokkos::ALL, 0);
  auto sqr_norm_j = Kokkos::subview(TMPNormView, Kokkos::ALL, 1);
  auto alpha      = Kokkos::subview(TMPNormView, Kokkos::ALL, 2);
  auto mask       = Kokkos::subview(TMPNormView, Kokkos::ALL, 3);
  auto tmp        = Kokkos::subview(TMPNormView, Kokkos::ALL, 4);

  TeamCopy<MemberType>::invoke(member, X, X_);
  // Deep copy of b into r_0:
  TeamCopy<MemberType>::invoke(member, B, R);

  // r_0 := b - A x_0
  member.team_barrier();
  A.template apply<Trans::NoTranspose, Mode::Team>(member, X_, R, -1, 1);
  member.team_barrier();

  // Deep copy of r_0 into p_0:
  TeamCopy<MemberType>::invoke(member, R, P);

  TeamDot<MemberType>::invoke(member, R, R, sqr_norm_0);
  member.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, numMatrices),
                       [&](const OrdinalType& i) { mask(i) = sqr_norm_0(i) > tolerance * tolerance ? 1. : 0; });

  TeamCopy1D::invoke(member, sqr_norm_0, sqr_norm_j);

  int status               = 1;
  int number_not_converged = 0;

  for (size_t j = 0; j < maximum_iteration; ++j) {
    // q := A p_j
    A.template apply<Trans::NoTranspose, Mode::Team>(member, P, Q);
    member.team_barrier();

    TeamDot<MemberType>::invoke(member, P, Q, tmp);
    member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, numMatrices),
                         [&](const OrdinalType& i) { alpha(i) = mask(i) != 0. ? sqr_norm_j(i) / tmp(i) : 0.; });
    member.team_barrier();

    // x_{j+1} := alpha p_j + x_j
    TeamAxpy<MemberType>::invoke(member, alpha, P, X_);
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
                         [&](const OrdinalType& i) { alpha(i) = mask(i) != 0. ? tmp(i) / sqr_norm_j(i) : 0.; });

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

  TeamCopy<MemberType>::invoke(member, X_, X);
  return status;
}

template <typename MemberType>
template <typename OperatorType, typename VectorViewType, typename KrylovHandleType>
KOKKOS_INLINE_FUNCTION int TeamCG<MemberType>::invoke(const MemberType& member, const OperatorType& A,
                                                      const VectorViewType& B, const VectorViewType& X,
                                                      const KrylovHandleType& handle) {
  const int strategy = handle.get_memory_strategy();
  if (strategy == 0) {
    using ScratchPadVectorViewType =
        Kokkos::View<typename VectorViewType::non_const_value_type**, typename VectorViewType::array_layout,
                     typename VectorViewType::execution_space::scratch_memory_space>;
    using ScratchPadNormViewType =
        Kokkos::View<typename KokkosKernels::ArithTraits<typename VectorViewType::non_const_value_type>::mag_type**,
                     typename VectorViewType::execution_space::scratch_memory_space>;

    const int numMatrices = X.extent(0);
    const int numRows     = X.extent(1);

    ScratchPadVectorViewType TMPView(member.team_scratch(handle.get_scratch_pad_level()), numMatrices, 4 * numRows);

    ScratchPadNormViewType TMPNormView(member.team_scratch(handle.get_scratch_pad_level()), numMatrices, 5);

    return invoke<OperatorType, VectorViewType, KrylovHandleType>(member, A, B, X, handle, TMPView, TMPNormView);
  }
  if (strategy == 1) {
    const int first_matrix = handle.first_index(member.league_rank());
    const int last_matrix  = handle.last_index(member.league_rank());

    using ScratchPadNormViewType =
        Kokkos::View<typename KokkosKernels::ArithTraits<typename VectorViewType::non_const_value_type>::mag_type**,
                     typename VectorViewType::execution_space::scratch_memory_space>;

    const int numMatrices = X.extent(0);

    auto TMPView = Kokkos::subview(handle.tmp_view, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);

    ScratchPadNormViewType TMPNormView(member.team_scratch(handle.get_scratch_pad_level()), numMatrices, 5);

    return invoke<OperatorType, VectorViewType, KrylovHandleType>(member, A, B, X, handle, TMPView, TMPNormView);
  }
  return 0;
}

}  // namespace KokkosBatched

#endif
