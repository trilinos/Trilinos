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

#ifndef KOKKOS_BLAS3_GEMM_DOTBASED_IMPL_HPP_
#define KOKKOS_BLAS3_GEMM_DOTBASED_IMPL_HPP_

#include "KokkosBlas_util.hpp"

namespace KokkosBlas {
namespace Impl {

// DotBasedGEMM implements the optimization for C = beta*C + alpha*A^TB
// with A and B matrices both being tall and skinny. C matrix is assumably
// small, so, each entry of C is computed by performing the dot product of
// respective columns of A and B matrices. Note that the dot products are
// performed on very long vectors, so, each dot product is distributed among
// numDivPerDot teams.

struct TagZero {};    // The init tag for beta=0
struct TagInit {};    // The init tag for beta!=0 and beta !=1
struct TagMult {};    // The multiplication tag for transposed A
struct TagMultCT {};  // The multiplication tag for conjugate-transposed A
template <class ExecSpace, class AV, class BV, class CV>
struct DotBasedGEMM {
  const AV A;
  const BV B;
  CV C;

  using scalar_A = typename AV::non_const_value_type;
  using size_A   = typename AV::size_type;
  using scalar_C = typename CV::non_const_value_type;
  using size_C   = typename CV::size_type;
  using AVT      = Kokkos::ArithTraits<scalar_A>;
  using CVT      = Kokkos::ArithTraits<scalar_C>;

  const scalar_A alpha;
  const scalar_C beta;

  const size_C numCrows;
  const size_C numCcols;

  size_C numDivPerDot;  // number of teams collectively performing a dot product
  size_C numTeams;      // total number of teams

  const size_A dotSize;  // the length of the vectors in the dot products

  DotBasedGEMM(const scalar_A& alpha_, const AV& A_, const BV& B_, const scalar_C& beta_, const CV& C_)
      : A(A_),
        B(B_),
        C(C_),
        alpha(alpha_),
        beta(beta_),
        numCrows(C.extent(0)),
        numCcols(C.extent(1)),
        dotSize(A.extent(0)) {}

  void run(const ExecSpace& space, bool conjugateTranspose) {
    multipleReductionWorkDistribution<ExecSpace, size_C>(dotSize, numCrows * numCcols, numDivPerDot);
    const size_C ndots = numCrows * numCcols;  // Number of dot products
    numTeams           = ndots * numDivPerDot;

    // Initialize C matrix if beta != 1
    if (beta == CVT::zero()) {
      Kokkos::MDRangePolicy<TagZero, ExecSpace, Kokkos::Rank<2>> policyInit(space, {0, 0}, {numCrows, numCcols});
      Kokkos::parallel_for("Initialize C for Dot Product Based GEMM", policyInit, *this);
    } else if (beta != CVT::one()) {
      Kokkos::MDRangePolicy<TagInit, ExecSpace, Kokkos::Rank<2>> policyInit(space, {0, 0}, {numCrows, numCcols});
      Kokkos::parallel_for("Initialize C for Dot Product Based GEMM", policyInit, *this);
    }

    // Multiply alpha*A^TB and add it to beta*C
    if (conjugateTranspose) {
      Kokkos::TeamPolicy<TagMultCT, ExecSpace> policyMult(space, numTeams, Kokkos::AUTO);
      Kokkos::parallel_for("Perform Dot Product Based GEMM", policyMult, *this);
    } else {
      Kokkos::TeamPolicy<TagMult, ExecSpace> policyMult(space, numTeams, Kokkos::AUTO);
      Kokkos::parallel_for("Perform Dot Product Based GEMM", policyMult, *this);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagZero&, const size_C& rowId, const size_C& colId) const { C(rowId, colId) = CVT::zero(); }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagInit&, const size_C& rowId, const size_C& colId) const {
    C(rowId, colId) = beta * C(rowId, colId);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagMult&, const typename Kokkos::TeamPolicy<ExecSpace>::member_type& teamMember) const {
    const size_C globalRank = teamMember.league_rank();
    const size_C localRank  = globalRank % numDivPerDot;
    const size_C i          = globalRank / numDivPerDot;
    const size_C rowId      = i / numCcols;
    const size_C colId      = i % numCcols;

    scalar_C result = CVT::zero();
    size_A begin    = localRank * (dotSize / numDivPerDot);
    size_A end      = (localRank + 1) * (dotSize / numDivPerDot);
    if (localRank == numDivPerDot - 1) end = dotSize;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(teamMember, begin, end),
        [&](const size_A k, scalar_C& update) { update += alpha * A(k, rowId) * B(k, colId); }, result);

    Kokkos::single(Kokkos::PerTeam(teamMember), [&]() { Kokkos::atomic_add(&C(rowId, colId), result); });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagMultCT&, const typename Kokkos::TeamPolicy<ExecSpace>::member_type& teamMember) const {
    const size_C globalRank = teamMember.league_rank();
    const size_C localRank  = globalRank % numDivPerDot;
    const size_C i          = globalRank / numDivPerDot;
    const size_C rowId      = i / numCcols;
    const size_C colId      = i % numCcols;

    scalar_C result = CVT::zero();
    size_A begin    = localRank * (dotSize / numDivPerDot);
    size_A end      = (localRank + 1) * (dotSize / numDivPerDot);
    if (localRank == numDivPerDot - 1) end = dotSize;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(teamMember, begin, end),
        [&](const size_A k, scalar_C& update) { update += alpha * AVT::conj(A(k, rowId)) * B(k, colId); }, result);

    Kokkos::single(Kokkos::PerTeam(teamMember), [&]() { Kokkos::atomic_add(&C(rowId, colId), result); });
  }
};

}  // namespace Impl
}  // namespace KokkosBlas

#endif
