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
#ifndef KOKKOS_BLAS2_MV_IMPL_GEMV_HPP_
#define KOKKOS_BLAS2_MV_IMPL_GEMV_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace KokkosBlas {
namespace Impl {

// Functor for a single-level parallel_for version of nontranspose
// GEMV.  The functor parallelizes over rows of the input matrix A.
template <class AViewType, class XViewType, class YViewType,
          const int alphaPreset,  // 0 or 1 are specific values; -1 means
                                  // general
          const int betaPreset,   // 0 or 1 are specific values; -1 means general
          class IndexType = typename AViewType::size_type>
struct SingleLevelNontransposeGEMV {
  using AlphaCoeffType = typename AViewType::non_const_value_type;
  using BetaCoeffType  = typename YViewType::non_const_value_type;
  using y_value_type   = typename YViewType::non_const_value_type;
  using AccumScalar    = typename std::conditional<std::is_same<y_value_type, Kokkos::Experimental::half_t>::value ||
                                                    std::is_same<y_value_type, Kokkos::Experimental::bhalf_t>::value,
                                                float, y_value_type>::type;

  SingleLevelNontransposeGEMV(const AlphaCoeffType& alpha, const AViewType& A, const XViewType& x,
                              const BetaCoeffType& beta, const YViewType& y)
      : alpha_(alpha), A_(A), x_(x), beta_(beta), y_(y) {
    static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<XViewType>::value, "XViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<YViewType>::value, "YViewType must be a Kokkos::View.");
    static_assert(static_cast<int>(AViewType::rank) == 2, "AViewType must have rank 2.");
    static_assert(static_cast<int>(XViewType::rank) == 1, "XViewType must have rank 1.");
    static_assert(static_cast<int>(YViewType::rank) == 1, "YViewType must have rank 1.");
    static_assert(std::is_integral<IndexType>::value, "IndexType must be an integer.");
    static_assert(alphaPreset == 0 || alphaPreset == 1 || alphaPreset == -1,
                  "Invalid alphaPreset value; valid values are 0, 1, and -1.");
    static_assert(betaPreset == 0 || betaPreset == 1 || betaPreset == -1,
                  "Invalid betaPreset value; valid values are 0, 1, and -1.");
  }

  // i is the current row of the input matrix A, and the current row
  // of the output vector y.
  KOKKOS_INLINE_FUNCTION void operator()(const IndexType& i) const {
    AccumScalar y_i;
    if (betaPreset == 0) {
      y_i = Kokkos::ArithTraits<AccumScalar>::zero();
    } else if (betaPreset == 1) {
      y_i = AccumScalar(y_(i));
    } else {  // beta_ != 0 and beta != 1
      y_i = beta_ * AccumScalar(y_(i));
    }

    const IndexType numCols = A_.extent(1);
    if (alphaPreset == 0) {
      ;  // do nothing
    } else if (alphaPreset == 1) {
      for (IndexType j = 0; j < numCols; ++j) {
        y_i += AccumScalar(A_(i, j)) * x_(j);
      }
    } else {  // alpha_ != 0 and alpha_ != 1
      for (IndexType j = 0; j < numCols; ++j) {
        y_i += alpha_ * AccumScalar(A_(i, j)) * x_(j);
      }
    }

    y_(i) = y_i;
  }

 private:
  AlphaCoeffType alpha_;
  typename AViewType::const_type A_;
  typename XViewType::const_type x_;
  BetaCoeffType beta_;
  YViewType y_;
};

// Functor for a single-level parallel_reduce version of (conjugate)
// transpose GEMV.  The functor parallelizes over rows of the input
// matrix A and the input vector x.  The output vector y is the
// reduction result.
//
// WARNING: NOT RECOMMENDED FOR GPU.  Reduction result may have
// arbitrary length.  This is bad on GPU because the GPU
// implementation of Kokkos::parallel_reduce may use shared memory for
// intermediate results.
template <class AViewType, class XViewType, class YViewType, const bool conj,
          const int alphaPreset,  // 0 or 1 are specific values; -1 means
                                  // general
          const int betaPreset,   // 0 or 1 are specific values; -1 means general
          class IndexType = typename AViewType::size_type>
struct SingleLevelTransposeGEMV {
  using y_value_type   = typename YViewType::non_const_value_type;
  using AlphaCoeffType = typename AViewType::non_const_value_type;
  using BetaCoeffType  = typename YViewType::non_const_value_type;
  using AccumScalar    = typename std::conditional<std::is_same<y_value_type, Kokkos::Experimental::half_t>::value ||
                                                    std::is_same<y_value_type, Kokkos::Experimental::bhalf_t>::value,
                                                float, y_value_type>::type;

  typedef AccumScalar value_type[];
  IndexType value_count;  // Kokkos needs this for reductions w/ array results

  SingleLevelTransposeGEMV(const AlphaCoeffType& alpha, const AViewType& A, const XViewType& x,
                           const BetaCoeffType& beta, const YViewType& y)
      : value_count(A.extent(1)), alpha_(alpha), A_(A), x_(x), beta_(beta), y_(y) {
    static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<XViewType>::value, "XViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<YViewType>::value, "YViewType must be a Kokkos::View.");
    static_assert(static_cast<int>(AViewType::rank) == 2, "AViewType must have rank 2.");
    static_assert(static_cast<int>(XViewType::rank) == 1, "XViewType must have rank 1.");
    static_assert(static_cast<int>(YViewType::rank) == 1, "YViewType must have rank 1.");
    static_assert(std::is_integral<IndexType>::value, "IndexType must be an integer.");
    static_assert(alphaPreset == 0 || alphaPreset == 1 || alphaPreset == -1,
                  "Invalid alphaPreset value; valid values are 0, 1, and -1.");
    static_assert(betaPreset == 0 || betaPreset == 1 || betaPreset == -1,
                  "Invalid betaPreset value; valid values are 0, 1, and -1.");
  }

 public:
  KOKKOS_INLINE_FUNCTION void init(value_type y_cur) const {
    for (IndexType j = 0; j < value_count; ++j) {
      y_cur[j] = Kokkos::ArithTraits<AccumScalar>::zero();
    }
  }

  KOKKOS_INLINE_FUNCTION void join(value_type dst, const value_type src) const {
    for (IndexType j = 0; j < value_count; ++j) {
      dst[j] += src[j];
    }
  }

  KOKKOS_INLINE_FUNCTION void final(value_type y_result) const {
    for (IndexType j = 0; j < value_count; ++j) {
      // Sum into initial y_ values; use beta as a pre-multiplier if nonzero.
      if (betaPreset == 0)
        y_(j) = y_value_type(alpha_ * y_result[j]);
      else
        y_(j) = y_value_type(beta_ * AccumScalar(y_(j)) + alpha_ * y_result[j]);
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(const IndexType& i, value_type y_cur) const {
    using Kokkos::ArithTraits;
    using KAT = ArithTraits<typename AViewType::non_const_value_type>;

    const auto x_i = x_(i);
    for (IndexType j = 0; j < value_count; ++j) {
      const auto A_ij = conj ? KAT::conj(A_(i, j)) : A_(i, j);
      y_cur[j] += AccumScalar(A_ij) * x_i;
    }
  }

 private:
  AlphaCoeffType alpha_;
  typename AViewType::const_type A_;
  typename XViewType::const_type x_;
  BetaCoeffType beta_;
  YViewType y_;
};

// Single-level parallel version of GEMV.
template <class ExecutionSpace, class AViewType, class XViewType, class YViewType,
          class IndexType = typename AViewType::size_type>
void singleLevelGemv(const ExecutionSpace& space, const char trans[], typename AViewType::const_value_type& alpha,
                     const AViewType& A, const XViewType& x, typename YViewType::const_value_type& beta,
                     const YViewType& y) {
  static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<XViewType>::value, "XViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<YViewType>::value, "YViewType must be a Kokkos::View.");
  static_assert(static_cast<int>(AViewType::rank) == 2, "AViewType must have rank 2.");
  static_assert(static_cast<int>(XViewType::rank) == 1, "XViewType must have rank 1.");
  static_assert(static_cast<int>(YViewType::rank) == 1, "YViewType must have rank 1.");
  static_assert(std::is_integral<IndexType>::value, "IndexType must be an integer");

  using y_value_type = typename YViewType::non_const_value_type;
  using policy_type  = Kokkos::RangePolicy<ExecutionSpace, IndexType>;

  using AlphaCoeffType = typename AViewType::non_const_value_type;
  using BetaCoeffType  = typename YViewType::non_const_value_type;

  policy_type range(space, 0, A.extent(0));
  const char tr = trans[0];

  // The transpose and conjugate transpose cases where A has zero rows
  // need special handling.  These are equivalent to y := beta*y.  We
  // could implement this using KokkosBlas::scal, but we don't want to
  // depend on that or its implementation details.  Instead, we reuse
  // an instantiation of the non-transpose case for alpha=0.
  if (A.extent(0) == 0 && (tr != 'N' && tr != 'n')) {
    if (beta == Kokkos::ArithTraits<BetaCoeffType>::zero()) {
      Kokkos::deep_copy(y, Kokkos::ArithTraits<BetaCoeffType>::zero());
    } else if (beta != Kokkos::ArithTraits<BetaCoeffType>::one()) {
      // "Fake out" a scal() by using the non-transpose alpha=0,
      // general beta case.  This assumes that the functor doesn't
      // check dimensions.
      using functor_type = SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType, 0, -1, IndexType>;
      functor_type functor(alpha, A, x, beta, y);
      Kokkos::parallel_for("KokkosBlas::gemv[SingleLevel]", policy_type(0, A.extent(1)), functor);
    }
    return;
  }

  if (tr == 'N' || tr == 'n') {
    if (alpha == Kokkos::ArithTraits<AlphaCoeffType>::zero()) {
      if (beta == Kokkos::ArithTraits<BetaCoeffType>::zero()) {
        // Fill y with zeros
        Kokkos::deep_copy(y, Kokkos::ArithTraits<y_value_type>::zero());
      } else if (beta == Kokkos::ArithTraits<BetaCoeffType>::one()) {
        // Do nothing (y := 1 * y)
      } else {  // beta != 0 && beta != 1
        using functor_type = SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType, 0, -1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_for("KokkosBlas::gemv[SingleLevel]", range, functor);
      }
    } else if (alpha == Kokkos::ArithTraits<AlphaCoeffType>::one()) {
      if (beta == Kokkos::ArithTraits<BetaCoeffType>::zero()) {
        using functor_type = SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType, 1, 0, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_for("KokkosBlas::gemv[SingleLevel]", range, functor);
      } else if (beta == Kokkos::ArithTraits<BetaCoeffType>::one()) {
        using functor_type = SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType, 1, 1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_for("KokkosBlas::gemv[SingleLevel]", range, functor);
      } else {  // beta != 0 && beta != 1
        using functor_type = SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType, 1, -1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_for("KokkosBlas::gemv[SingleLevel]", range, functor);
      }
    } else {  // alpha != 0 and alpha != 1
      if (beta == Kokkos::ArithTraits<BetaCoeffType>::zero()) {
        using functor_type = SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType, -1, 0, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_for("KokkosBlas::gemv[SingleLevel]", range, functor);
      } else if (beta == Kokkos::ArithTraits<BetaCoeffType>::one()) {
        using functor_type = SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType, -1, 1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_for("KokkosBlas::gemv[SingleLevel]", range, functor);
      } else {  // beta != 0 && beta != 1
        using functor_type = SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType, -1, -1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_for("KokkosBlas::gemv[SingleLevel]", range, functor);
      }
    }
  } else if (tr == 'T' || tr == 't') {  // transpose, no conjugate
    if (alpha == Kokkos::ArithTraits<AlphaCoeffType>::zero()) {
      if (beta == Kokkos::ArithTraits<BetaCoeffType>::zero()) {
        // Fill y with zeros
        Kokkos::deep_copy(y, Kokkos::ArithTraits<y_value_type>::zero());
      } else if (beta == Kokkos::ArithTraits<BetaCoeffType>::one()) {
        // Do nothing (y := 1 * y)
      } else {  // beta != 0 && beta != 1
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, false, 0, -1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      }
    } else if (alpha == Kokkos::ArithTraits<AlphaCoeffType>::one()) {
      if (beta == Kokkos::ArithTraits<BetaCoeffType>::zero()) {
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, false, 1, 0, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      } else if (beta == Kokkos::ArithTraits<BetaCoeffType>::one()) {
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, false, 1, 1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      } else {  // beta != 0 && beta != 1
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, false, 1, -1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      }
    } else {  // alpha != 0 and alpha != 1
      if (beta == Kokkos::ArithTraits<BetaCoeffType>::zero()) {
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, false, -1, 0, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      } else if (beta == Kokkos::ArithTraits<BetaCoeffType>::one()) {
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, false, -1, 1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      } else {  // beta != 0 && beta != 1
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, false, -1, -1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      }
    }
  } else if (tr == 'C' || tr == 'c' || tr == 'H' || tr == 'h') {  // conj xpose
    if (alpha == Kokkos::ArithTraits<AlphaCoeffType>::zero()) {
      if (beta == Kokkos::ArithTraits<BetaCoeffType>::zero()) {
        // Fill y with zeros
        Kokkos::deep_copy(y, Kokkos::ArithTraits<y_value_type>::zero());
      } else if (beta == Kokkos::ArithTraits<BetaCoeffType>::one()) {
        // Do nothing (y := 1 * y)
      } else {  // beta != 0 && beta != 1
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, true, 0, -1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      }
    } else if (alpha == Kokkos::ArithTraits<AlphaCoeffType>::one()) {
      if (beta == Kokkos::ArithTraits<BetaCoeffType>::zero()) {
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, true, 1, 0, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      } else if (beta == Kokkos::ArithTraits<BetaCoeffType>::one()) {
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, true, 1, 1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      } else {  // beta != 0 && beta != 1
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, true, 1, -1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      }
    } else {  // alpha != 0 and alpha != 1
      if (beta == Kokkos::ArithTraits<BetaCoeffType>::zero()) {
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, true, -1, 0, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      } else if (beta == Kokkos::ArithTraits<BetaCoeffType>::one()) {
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, true, -1, 1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      } else {  // beta != 0 && beta != 1
        using functor_type = SingleLevelTransposeGEMV<AViewType, XViewType, YViewType, true, -1, -1, IndexType>;
        functor_type functor(alpha, A, x, beta, y);
        Kokkos::parallel_reduce("KokkosBlas::gemv[SingleLevelTranspose]", range, functor);
      }
    }
  }
}

struct TwoLevelGEMV_LayoutLeftTag {};
struct TwoLevelGEMV_LayoutRightTag {};

// ---------------------------------------------------------------------------------------------
// Functor for a two-level parallel_reduce version of GEMV (non-transpose),
// designed for performance on GPU. Kernel depends on the layout of A.
template <class ExecutionSpace, class AViewType, class XViewType, class YViewType,
          class IndexType = typename AViewType::size_type>
struct TwoLevelGEMV {
  using y_value_type   = typename YViewType::non_const_value_type;
  using AlphaCoeffType = typename AViewType::non_const_value_type;
  using BetaCoeffType  = typename YViewType::non_const_value_type;
  using AccumScalar    = typename std::conditional<std::is_same<y_value_type, Kokkos::Experimental::half_t>::value ||
                                                    std::is_same<y_value_type, Kokkos::Experimental::bhalf_t>::value,
                                                float, y_value_type>::type;

  using policy_type = Kokkos::TeamPolicy<ExecutionSpace>;
  using member_type = typename policy_type::member_type;

  TwoLevelGEMV(const AlphaCoeffType& alpha, const AViewType& A, const XViewType& x, const BetaCoeffType& beta,
               const YViewType& y)
      : alpha_(alpha), A_(A), x_(x), beta_(beta), y_(y) {
    static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<XViewType>::value, "XViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<YViewType>::value, "YViewType must be a Kokkos::View.");
    static_assert(static_cast<int>(AViewType::rank) == 2, "AViewType must have rank 2.");
    static_assert(static_cast<int>(XViewType::rank) == 1, "XViewType must have rank 1.");
    static_assert(static_cast<int>(YViewType::rank) == 1, "YViewType must have rank 1.");
    static_assert(std::is_integral<IndexType>::value, "IndexType must be an integer.");
  }

 public:
  // LayoutLeft version: 32xK blocks.
  //  -Each team handles block rows.
  //  -Groups of 32 threads handle N/teamsize columns sequentially, placing
  //  results into shared. -Then individual thread results are combined with
  //  parallel_reduce.
  KOKKOS_INLINE_FUNCTION void operator()(TwoLevelGEMV_LayoutLeftTag, const member_type& team) const {
    using KAT  = Kokkos::ArithTraits<y_value_type>;
    using AKAT = Kokkos::ArithTraits<AccumScalar>;
    // Allocate a Scalar in shared for each thread
    AccumScalar* blockResult = (AccumScalar*)team.team_shmem().get_shmem(32 * sizeof(AccumScalar));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 32), [&](int i) { blockResult[i] = AKAT::zero(); });
    team.team_barrier();
    // Which block this thread will work on
    int block = team.team_rank() / 32;
    // Which row in the block this thread will work on
    IndexType row           = team.league_rank() * 32 + team.team_rank() % 32;
    IndexType blockColStart = columnsPerThread * block;
    AccumScalar localSum    = AKAT::zero();
    // compute local sum
    if (row < (IndexType)A_.extent(0)) {
      for (IndexType col = blockColStart; col < blockColStart + columnsPerThread && col < A_.extent(1); col++) {
        // A access is coalesced, x access is a broadcast
        localSum += AccumScalar(A_(row, col)) * AccumScalar(x_(col));
      }
      // atomically combine local result into shared
      Kokkos::atomic_add(&blockResult[team.team_rank() % 32], localSum);
    }
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 32), [&](int i) {
      IndexType yrow = team.league_rank() * 32 + i;
      if (yrow < (IndexType)A_.extent(0)) {
        if (beta_ == KAT::zero())
          y_(yrow) = y_value_type(alpha_ * blockResult[i]);
        else
          y_(yrow) = y_value_type(beta_ * AccumScalar(y_(yrow)) + alpha_ * blockResult[i]);
      }
    });
  }

  // LayoutRight version: one team per row
  KOKKOS_INLINE_FUNCTION void operator()(TwoLevelGEMV_LayoutRightTag, const member_type& team) const {
    using KAT = Kokkos::ArithTraits<y_value_type>;

    const IndexType N = A_.extent(1);
    const int i       = team.league_rank();  // batch id

    // parallel-reduce to compute val += A(:,j)' * x
    AccumScalar val;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, N),
        [&](const int j, AccumScalar& update) { update += AccumScalar(A_(i, j)) * x_(j); }, val);

    // compute yj = beta*yj + alpha*val
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      if (beta_ == KAT::zero())
        y_(i) = y_value_type(alpha_ * val);
      else
        y_(i) = y_value_type(beta_ * AccumScalar(y_(i)) + alpha_ * val);
    });
  }

  IndexType columnsPerThread;

 private:
  AlphaCoeffType alpha_;
  typename AViewType::const_type A_;
  typename XViewType::const_type x_;
  BetaCoeffType beta_;
  YViewType y_;
};

// ---------------------------------------------------------------------------------------------
// Functor for a two-level parallel_reduce version of (conjugate)
// transpose GEMV.  The functor uses parallel-for over the columns of the input
// matrix A and each team uses parallel-reduce over the row of its column.
// The output vector y is the reduction result.
template <class ExecutionSpace, class AViewType, class XViewType, class YViewType, const bool conj,
          class IndexType = typename AViewType::size_type>
struct TwoLevelTransposeGEMV {
  using y_value_type   = typename YViewType::non_const_value_type;
  using AlphaCoeffType = typename AViewType::non_const_value_type;
  using BetaCoeffType  = typename YViewType::non_const_value_type;
  using AccumScalar    = typename std::conditional<std::is_same<y_value_type, Kokkos::Experimental::half_t>::value ||
                                                    std::is_same<y_value_type, Kokkos::Experimental::bhalf_t>::value,
                                                float, y_value_type>::type;

  using policy_type = Kokkos::TeamPolicy<ExecutionSpace>;
  using member_type = typename policy_type::member_type;

  TwoLevelTransposeGEMV(const AlphaCoeffType& alpha, const AViewType& A, const XViewType& x, const BetaCoeffType& beta,
                        const YViewType& y)
      : alpha_(alpha), A_(A), x_(x), beta_(beta), y_(y) {
    static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<XViewType>::value, "XViewType must be a Kokkos::View.");
    static_assert(Kokkos::is_view<YViewType>::value, "YViewType must be a Kokkos::View.");
    static_assert(static_cast<int>(AViewType::rank) == 2, "AViewType must have rank 2.");
    static_assert(static_cast<int>(XViewType::rank) == 1, "XViewType must have rank 1.");
    static_assert(static_cast<int>(YViewType::rank) == 1, "YViewType must have rank 1.");
    static_assert(std::is_integral<IndexType>::value, "IndexType must be an integer.");
  }

 public:
  KOKKOS_INLINE_FUNCTION void operator()(const member_type& team) const {
    using Kokkos::ArithTraits;
    using KAT_A = ArithTraits<typename AViewType::non_const_value_type>;
    using KAT_Y = ArithTraits<typename YViewType::non_const_value_type>;

    const IndexType M = A_.extent(0);
    const int j       = team.league_rank();  // batch id

    // parallel-reduce to compute val += A(:,j)' * x
    AccumScalar val;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, M),
        [&](const int i, AccumScalar& update) {
          const auto x_i  = x_(i);
          const auto A_ij = conj ? KAT_A::conj(A_(i, j)) : A_(i, j);
          update += AccumScalar(A_ij) * x_i;
        },
        val);

    // compute yj = beta*yj + alpha*val
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      if (beta_ == KAT_Y::zero())
        y_(j) = y_value_type(alpha_ * val);
      else
        y_(j) = y_value_type(beta_ * AccumScalar(y_(j)) + alpha_ * val);
    });
  }

 private:
  AlphaCoeffType alpha_;
  typename AViewType::const_type A_;
  typename XViewType::const_type x_;
  BetaCoeffType beta_;
  YViewType y_;
};

// Two-level parallel version of GEMV.
template <class ExecutionSpace, class AViewType, class XViewType, class YViewType,
          class IndexType = typename AViewType::size_type>
void twoLevelGemv(const ExecutionSpace& space, const char trans[], typename AViewType::const_value_type& alpha,
                  const AViewType& A, const XViewType& x, typename YViewType::const_value_type& beta,
                  const YViewType& y) {
  static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<XViewType>::value, "XViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<YViewType>::value, "YViewType must be a Kokkos::View.");
  static_assert(static_cast<int>(AViewType::rank) == 2, "AViewType must have rank 2.");
  static_assert(static_cast<int>(XViewType::rank) == 1, "XViewType must have rank 1.");
  static_assert(static_cast<int>(YViewType::rank) == 1, "YViewType must have rank 1.");
  static_assert(std::is_integral<IndexType>::value, "IndexType must be an integer");

  using y_value_type      = typename YViewType::non_const_value_type;
  using team_policy_type  = Kokkos::TeamPolicy<ExecutionSpace>;
  using range_policy_type = Kokkos::RangePolicy<ExecutionSpace, IndexType>;

  using Kokkos::ArithTraits;
  using KAT  = ArithTraits<typename AViewType::non_const_value_type>;
  using YKAT = ArithTraits<typename YViewType::non_const_value_type>;

  const char tr = toupper(trans[0]);

  // The transpose and conjugate transpose cases where A has zero rows
  // need special handling.  These are equivalent to y := beta*y.  We
  // could implement this using KokkosBlas::scal, but we don't want to
  // depend on that or its implementation details.  Instead, we reuse
  // an instantiation of the non-transpose case for alpha=0.
  if (y.extent(0) == 0) {
    // no entries to update
    return;
  } else if (x.extent(0) == 0) {
    if (beta == YKAT::zero()) {
      Kokkos::deep_copy(y, KAT::zero());
    } else if (beta != YKAT::one()) {
      // "Fake out" a scal() by using the non-transpose alpha=0,
      // general beta case.  This assumes that the functor doesn't
      // check dimensions.
      using functor_type = SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType, 0, -1, IndexType>;
      functor_type functor(alpha, A, x, beta, y);
      Kokkos::parallel_for("KokkosBlas::gemv[SingleLevel]", range_policy_type(space, 0, y.extent(0)), functor);
    }
    return;
  }

  if (tr == 'N') {
    constexpr bool isLayoutLeft = std::is_same<typename AViewType::array_layout, Kokkos::LayoutLeft>::value;
    // Both kernels work for both layouts - the only difference is access
    // pattern.
    using layout_tag =
        typename std::conditional<isLayoutLeft, TwoLevelGEMV_LayoutLeftTag, TwoLevelGEMV_LayoutRightTag>::type;
    using tagged_policy = Kokkos::TeamPolicy<ExecutionSpace, layout_tag>;
    using functor_type  = TwoLevelGEMV<ExecutionSpace, AViewType, XViewType, YViewType, IndexType>;
    functor_type functor(alpha, A, x, beta, y);
    tagged_policy team;
    if constexpr (isLayoutLeft) {
      using AccumScalar =
          typename std::conditional<std::is_same<y_value_type, Kokkos::Experimental::half_t>::value ||
                                        std::is_same<y_value_type, Kokkos::Experimental::bhalf_t>::value,
                                    float, y_value_type>::type;
      size_t sharedPerTeam = 32 * sizeof(AccumScalar);
      IndexType numTeams   = (A.extent(0) + 31) / 32;
      tagged_policy temp(space, 1, 1);
      temp.set_scratch_size(0, Kokkos::PerTeam(sharedPerTeam));
      int teamSize = temp.team_size_recommended(functor, Kokkos::ParallelForTag());
      // make sure teamSize is a multiple of 32
      teamSize -= teamSize % 32;
      // don't make teamSize larger than what's useful
      if ((size_t)teamSize > 32 * A.extent(1)) teamSize = 32 * A.extent(1);
        // FIXME SYCL: team_size_recommended() returns too big of a team size.
        // Kernel hangs with 1024 threads on XEHP.
#ifdef KOKKOS_ENABLE_SYCL
      if (std::is_same<ExecutionSpace, Kokkos::Experimental::SYCL>::value) {
        if (teamSize > 256) teamSize = 256;
      }
#endif
      int numBlocks            = teamSize / 32;
      functor.columnsPerThread = (A.extent(1) + numBlocks - 1) / numBlocks;
      team = tagged_policy(space, numTeams, teamSize).set_scratch_size(0, Kokkos::PerTeam(sharedPerTeam));
    } else {
      // LayoutRight: one team per row
      team = tagged_policy(space, A.extent(0), Kokkos::AUTO);
    }
    Kokkos::parallel_for("KokkosBlas::gemv[twoLevel]", team, functor);
  } else {
    if (alpha == KAT::zero() && beta == KAT::zero()) {
      // Fill y with zeros
      Kokkos::deep_copy(y, Kokkos::ArithTraits<y_value_type>::zero());
    } else if (alpha == KAT::zero() && beta == KAT::one()) {
      // Do nothing (y := 1 * y)
    } else if (tr == 'T') {
      // transpose, and not conj transpose
      team_policy_type team(space, A.extent(1), Kokkos::AUTO);
      using functor_type = TwoLevelTransposeGEMV<ExecutionSpace, AViewType, XViewType, YViewType, false, IndexType>;
      functor_type functor(alpha, A, x, beta, y);
      Kokkos::parallel_for("KokkosBlas::gemv[twoLevelTranspose]", team, functor);
    } else if (tr == 'C' || tr == 'H') {
      // conjugate transpose
      team_policy_type team(space, A.extent(1), Kokkos::AUTO);
      using functor_type = TwoLevelTransposeGEMV<ExecutionSpace, AViewType, XViewType, YViewType, true, IndexType>;
      functor_type functor(alpha, A, x, beta, y);
      Kokkos::parallel_for("KokkosBlas::gemv[twoLevelTranspose]", team, functor);
    }
  }
}

// generalGemv: use 1 level (Range) or 2 level (Team) implementation,
// depending on whether execution space is CPU or GPU. enable_if makes sure
// unused kernels are not instantiated.
template <class ExecutionSpace, class AViewType, class XViewType, class YViewType, class IndexType,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>()>::type* = nullptr>
void generalGemvImpl(const ExecutionSpace& space, const char trans[], typename AViewType::const_value_type& alpha,
                     const AViewType& A, const XViewType& x, typename YViewType::const_value_type& beta,
                     const YViewType& y) {
  singleLevelGemv(space, trans, alpha, A, x, beta, y);
}

template <class ExecutionSpace, class AViewType, class XViewType, class YViewType, class IndexType,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>()>::type* = nullptr>
void generalGemvImpl(const ExecutionSpace& space, const char trans[], typename AViewType::const_value_type& alpha,
                     const AViewType& A, const XViewType& x, typename YViewType::const_value_type& beta,
                     const YViewType& y) {
  twoLevelGemv(space, trans, alpha, A, x, beta, y);
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOS_BLAS2_MV_IMPL_GEMV_HPP_
