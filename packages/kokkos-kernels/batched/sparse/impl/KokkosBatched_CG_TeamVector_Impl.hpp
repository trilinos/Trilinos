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
#ifndef __KOKKOSBATCHED_CG_TEAMVECTOR_IMPL_HPP__
#define __KOKKOSBATCHED_CG_TEAMVECTOR_IMPL_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Axpy.hpp"
#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBatched_Dot.hpp"
#include "KokkosBatched_Spmv.hpp"
#include "KokkosBatched_Xpay.hpp"

namespace KokkosBatched {

///
/// TeamVector CG
///   Two nested parallel_for with both TeamVectorRange and ThreadVectorRange
///   (or one with TeamVectorRange) are used inside.
///

template <typename MemberType>
template <typename OperatorType, typename VectorViewType, typename KrylovHandleType, typename TMPViewType,
          typename TMPNormViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorCG<MemberType>::invoke(const MemberType& member, const OperatorType& A,
                                                            const VectorViewType& _B, const VectorViewType& _X,
                                                            const KrylovHandleType& handle, const TMPViewType& _TMPView,
                                                            const TMPNormViewType& _TMPNormView) {
  typedef int OrdinalType;
  typedef typename Kokkos::ArithTraits<typename VectorViewType::non_const_value_type>::mag_type MagnitudeType;

  const size_t maximum_iteration = handle.get_max_iteration();
  const MagnitudeType tolerance  = handle.get_tolerance();

  using TeamVectorCopy1D = TeamVectorCopy<MemberType, Trans::NoTranspose, 1>;

  const OrdinalType numMatrices = _X.extent(0);
  const OrdinalType numRows     = _X.extent(1);

  int offset_P = 0;
  int offset_Q = offset_P + numRows;
  int offset_R = offset_Q + numRows;
  int offset_X = offset_R + numRows;

  auto P = Kokkos::subview(_TMPView, Kokkos::ALL, Kokkos::make_pair(offset_P, offset_P + numRows));
  auto Q = Kokkos::subview(_TMPView, Kokkos::ALL, Kokkos::make_pair(offset_Q, offset_Q + numRows));
  auto R = Kokkos::subview(_TMPView, Kokkos::ALL, Kokkos::make_pair(offset_R, offset_R + numRows));
  auto X = Kokkos::subview(_TMPView, Kokkos::ALL, Kokkos::make_pair(offset_X, offset_X + numRows));

  auto sqr_norm_0 = Kokkos::subview(_TMPNormView, Kokkos::ALL, 0);
  auto sqr_norm_j = Kokkos::subview(_TMPNormView, Kokkos::ALL, 1);
  auto alpha      = Kokkos::subview(_TMPNormView, Kokkos::ALL, 2);
  auto mask       = Kokkos::subview(_TMPNormView, Kokkos::ALL, 3);
  auto tmp        = Kokkos::subview(_TMPNormView, Kokkos::ALL, 4);

  TeamVectorCopy<MemberType>::invoke(member, _X, X);
  // Deep copy of b into r_0:
  TeamVectorCopy<MemberType>::invoke(member, _B, R);

  // r_0 := b - A x_0
  member.team_barrier();
  A.template apply<Trans::NoTranspose, Mode::TeamVector>(member, X, R, -1, 1);
  member.team_barrier();

  // Deep copy of r_0 into p_0:
  TeamVectorCopy<MemberType>::invoke(member, R, P);

  TeamVectorDot<MemberType>::invoke(member, R, R, sqr_norm_0);
  member.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, numMatrices),
                       [&](const OrdinalType& i) { mask(i) = sqr_norm_0(i) > tolerance * tolerance ? 1. : 0; });

  TeamVectorCopy1D::invoke(member, sqr_norm_0, sqr_norm_j);

  int status               = 1;
  int number_not_converged = 0;

  for (size_t j = 0; j < maximum_iteration; ++j) {
    // q := A p_j
    A.template apply<Trans::NoTranspose, Mode::TeamVector>(member, P, Q);
    member.team_barrier();

    TeamVectorDot<MemberType>::invoke(member, P, Q, tmp);
    member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, numMatrices),
                         [&](const OrdinalType& i) { alpha(i) = mask(i) != 0. ? sqr_norm_j(i) / tmp(i) : 0.; });
    member.team_barrier();

    // x_{j+1} := alpha p_j + x_j
    TeamVectorAxpy<MemberType>::invoke(member, alpha, P, X);
    member.team_barrier();

    // r_{j+1} := - alpha q + r_j
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, numMatrices),
                         [&](const OrdinalType& i) { alpha(i) = -alpha(i); });
    member.team_barrier();

    TeamVectorAxpy<MemberType>::invoke(member, alpha, Q, R);
    member.team_barrier();

    TeamVectorDot<MemberType>::invoke(member, R, R, tmp);
    member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, numMatrices),
                         [&](const OrdinalType& i) { alpha(i) = mask(i) != 0. ? tmp(i) / sqr_norm_j(i) : 0.; });

    TeamVectorCopy1D::invoke(member, tmp, sqr_norm_j);

    // Relative convergence check:
    number_not_converged = 0;
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(member, 0, numMatrices),
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
    TeamVectorXpay<MemberType>::invoke(member, alpha, R, P);
    member.team_barrier();
  }

  TeamVectorCopy<MemberType>::invoke(member, X, _X);
  return status;
}

template <typename MemberType>
template <typename OperatorType, typename VectorViewType, typename KrylovHandleType>
KOKKOS_INLINE_FUNCTION int TeamVectorCG<MemberType>::invoke(const MemberType& member, const OperatorType& A,
                                                            const VectorViewType& _B, const VectorViewType& _X,
                                                            const KrylovHandleType& handle) {
  const int strategy = handle.get_memory_strategy();
  if (strategy == 0) {
    using ScratchPadVectorViewType =
        Kokkos::View<typename VectorViewType::non_const_value_type**, typename VectorViewType::array_layout,
                     typename VectorViewType::execution_space::scratch_memory_space>;
    using ScratchPadNormViewType =
        Kokkos::View<typename Kokkos::ArithTraits<typename VectorViewType::non_const_value_type>::mag_type**,
                     typename VectorViewType::execution_space::scratch_memory_space>;

    const int numMatrices = _X.extent(0);
    const int numRows     = _X.extent(1);

    ScratchPadVectorViewType _TMPView(member.team_scratch(handle.get_scratch_pad_level()), numMatrices, 4 * numRows);

    ScratchPadNormViewType _TMPNormView(member.team_scratch(handle.get_scratch_pad_level()), numMatrices, 5);

    return invoke<OperatorType, VectorViewType, KrylovHandleType>(member, A, _B, _X, handle, _TMPView, _TMPNormView);
  }
  if (strategy == 1) {
    const int first_matrix = handle.first_index(member.league_rank());
    const int last_matrix  = handle.last_index(member.league_rank());

    using ScratchPadNormViewType =
        Kokkos::View<typename Kokkos::ArithTraits<typename VectorViewType::non_const_value_type>::mag_type**,
                     typename VectorViewType::execution_space::scratch_memory_space>;

    const int numMatrices = _X.extent(0);

    auto _TMPView = Kokkos::subview(handle.tmp_view, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);

    ScratchPadNormViewType _TMPNormView(member.team_scratch(handle.get_scratch_pad_level()), numMatrices, 5);

    return invoke<OperatorType, VectorViewType, KrylovHandleType>(member, A, _B, _X, handle, _TMPView, _TMPNormView);
  }
  return 0;
}

}  // namespace KokkosBatched

#endif
