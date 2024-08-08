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

#ifndef KOKKOSBLAS2_SYR2_IMPL_HPP_
#define KOKKOSBLAS2_SYR2_IMPL_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace KokkosBlas {
namespace Impl {

// Functor for the thread parallel version of SYR2.
// This functor parallelizes over rows of the input matrix A.
template <class XViewType, class YViewType, class AViewType, class IndexType, bool tJustTranspose, bool tJustUp>
struct ThreadParallelSYR2 {
  using AlphaCoeffType = typename AViewType::non_const_value_type;
  using XComponentType = typename XViewType::non_const_value_type;
  using YComponentType = typename YViewType::non_const_value_type;
  using AComponentType = typename AViewType::non_const_value_type;

  ThreadParallelSYR2(const AlphaCoeffType& alpha, const XViewType& x, const YViewType& y, const AViewType& A)
      : alpha_(alpha), x_(x), y_(y), A_(A) {
    // Nothing to do
  }

  KOKKOS_INLINE_FUNCTION void operator()(const IndexType& i) const {
    if (alpha_ == Kokkos::ArithTraits<AlphaCoeffType>::zero()) {
      // Nothing to do
    } else if ((x_(i) == Kokkos::ArithTraits<XComponentType>::zero()) &&
               (y_(i) == Kokkos::ArithTraits<YComponentType>::zero())) {
      // Nothing to do
    } else {
      const XComponentType x_fixed(x_(i));
      const YComponentType y_fixed(y_(i));
      const IndexType N(A_.extent(1));

      if constexpr (tJustTranspose) {
        if (x_fixed != Kokkos::ArithTraits<XComponentType>::zero()) {
          for (IndexType j = 0; j < N; ++j) {
            if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
              A_(i, j) += AComponentType(alpha_ * x_fixed * y_(j));
            }
          }
        }
        if (y_fixed != Kokkos::ArithTraits<YComponentType>::zero()) {
          for (IndexType j = 0; j < N; ++j) {
            if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
              A_(i, j) += AComponentType(alpha_ * y_fixed * x_(j));
            }
          }
        }
      } else {
        if (x_fixed != Kokkos::ArithTraits<XComponentType>::zero()) {
          for (IndexType j = 0; j < N; ++j) {
            if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
              A_(i, j) += AComponentType(alpha_ * x_fixed * Kokkos::ArithTraits<YComponentType>::conj(y_(j)));
            }
          }
        }
        if (y_fixed != Kokkos::ArithTraits<YComponentType>::zero()) {
          for (IndexType j = 0; j < N; ++j) {
            if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
              A_(i, j) += AComponentType(Kokkos::ArithTraits<AlphaCoeffType>::conj(alpha_) * y_fixed *
                                         Kokkos::ArithTraits<XComponentType>::conj(x_(j)));
            }
          }
        }
      }
    }
  }

 private:
  AlphaCoeffType alpha_;
  typename XViewType::const_type x_;
  typename YViewType::const_type y_;
  AViewType A_;
};

// Thread parallel version of SYR2.
template <class ExecutionSpace, class XViewType, class YViewType, class AViewType, class IndexType, bool tJustTranspose,
          bool tJustUp>
void threadParallelSyr2(const ExecutionSpace& space, const typename AViewType::const_value_type& alpha,
                        const XViewType& x, const YViewType& y, const AViewType& A) {
  static_assert(std::is_integral<IndexType>::value, "IndexType must be an integer");

  using AlphaCoeffType = typename AViewType::non_const_value_type;

  if (x.extent(0) == 0) {
    // no entries to update
  } else if (y.extent(0) == 0) {
    // no entries to update
  } else if (alpha == Kokkos::ArithTraits<AlphaCoeffType>::zero()) {
    // no entries to update
  } else {
    Kokkos::RangePolicy<ExecutionSpace, IndexType> rangePolicy(space, 0, A.extent(0));
    ThreadParallelSYR2<XViewType, YViewType, AViewType, IndexType, tJustTranspose, tJustUp> functor(alpha, x, y, A);
    Kokkos::parallel_for("KokkosBlas::syr2[threadParallel]", rangePolicy, functor);
  }
}

struct TeamParallelSYR2_LayoutLeftTag {};
struct TeamParallelSYR2_LayoutRightTag {};

// ---------------------------------------------------------------------------------------------

// Functor for the team parallel version of SYR2, designed for
// performance on GPUs. The kernel depends on the layout of A.
template <class ExecutionSpace, class XViewType, class YViewType, class AViewType, class IndexType, bool tJustTranspose,
          bool tJustUp>
struct TeamParallelSYR2 {
  using AlphaCoeffType = typename AViewType::non_const_value_type;
  using XComponentType = typename XViewType::non_const_value_type;
  using YComponentType = typename YViewType::non_const_value_type;
  using AComponentType = typename AViewType::non_const_value_type;

  using policy_type = Kokkos::TeamPolicy<ExecutionSpace>;
  using member_type = typename policy_type::member_type;

  TeamParallelSYR2(const AlphaCoeffType& alpha, const XViewType& x, const YViewType& y, const AViewType& A)
      : alpha_(alpha), x_(x), y_(y), A_(A) {
    // Nothing to do
  }

 public:
  // LayoutLeft version: one team per column
  KOKKOS_INLINE_FUNCTION void operator()(TeamParallelSYR2_LayoutLeftTag, const member_type& team) const {
    if (alpha_ == Kokkos::ArithTraits<AlphaCoeffType>::zero()) {
      // Nothing to do
    } else {
      const IndexType j(team.league_rank());
      if ((x_(j) == Kokkos::ArithTraits<XComponentType>::zero()) &&
          (y_(j) == Kokkos::ArithTraits<YComponentType>::zero())) {
        // Nothing to do
      } else {
        const IndexType M(A_.extent(0));
        if constexpr (tJustTranspose) {
          const XComponentType x_fixed(x_(j));
          const YComponentType y_fixed(y_(j));
          if (y_fixed != Kokkos::ArithTraits<YComponentType>::zero()) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, M), [&](const IndexType& i) {
              if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
                A_(i, j) += AComponentType(alpha_ * x_(i) * y_fixed);
              }
            });
          }
          if (x_fixed != Kokkos::ArithTraits<XComponentType>::zero()) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, M), [&](const IndexType& i) {
              if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
                A_(i, j) += AComponentType(alpha_ * y_(i) * x_fixed);
              }
            });
          }
        } else {
          const XComponentType x_fixed(Kokkos::ArithTraits<XComponentType>::conj(x_(j)));
          const YComponentType y_fixed(Kokkos::ArithTraits<YComponentType>::conj(y_(j)));
          if (y_fixed != Kokkos::ArithTraits<YComponentType>::zero()) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, M), [&](const IndexType& i) {
              if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
                A_(i, j) += AComponentType(alpha_ * x_(i) * y_fixed);
              }
            });
          }
          if (x_fixed != Kokkos::ArithTraits<XComponentType>::zero()) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, M), [&](const IndexType& i) {
              if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
                A_(i, j) += AComponentType(Kokkos::ArithTraits<AlphaCoeffType>::conj(alpha_) * y_(i) * x_fixed);
              }
            });
          }
        }
      }
    }
  }

  // LayoutRight version: one team per row
  KOKKOS_INLINE_FUNCTION void operator()(TeamParallelSYR2_LayoutRightTag, const member_type& team) const {
    if (alpha_ == Kokkos::ArithTraits<AlphaCoeffType>::zero()) {
      // Nothing to do
    } else {
      const IndexType i(team.league_rank());
      if ((x_(i) == Kokkos::ArithTraits<XComponentType>::zero()) &&
          (y_(i) == Kokkos::ArithTraits<YComponentType>::zero())) {
        // Nothing to do
      } else {
        const IndexType N(A_.extent(1));
        const XComponentType x_fixed(x_(i));
        const YComponentType y_fixed(y_(i));
        if constexpr (tJustTranspose) {
          if (x_fixed != Kokkos::ArithTraits<XComponentType>::zero()) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const IndexType& j) {
              if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
                A_(i, j) += AComponentType(alpha_ * x_fixed * y_(j));
              }
            });
          }
          if (y_fixed != Kokkos::ArithTraits<YComponentType>::zero()) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const IndexType& j) {
              if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
                A_(i, j) += AComponentType(alpha_ * y_fixed * x_(j));
              }
            });
          }
        } else {
          if (x_fixed != Kokkos::ArithTraits<XComponentType>::zero()) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const IndexType& j) {
              if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
                A_(i, j) += AComponentType(alpha_ * x_fixed * Kokkos::ArithTraits<YComponentType>::conj(y_(j)));
              }
            });
          }
          if (y_fixed != Kokkos::ArithTraits<YComponentType>::zero()) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const IndexType& j) {
              if (((tJustUp == true) && (i <= j)) || ((tJustUp == false) && (i >= j))) {
                A_(i, j) += AComponentType(Kokkos::ArithTraits<AlphaCoeffType>::conj(alpha_) * y_fixed *
                                           Kokkos::ArithTraits<XComponentType>::conj(x_(j)));
              }
            });
          }
        }
      }
    }
  }

 private:
  AlphaCoeffType alpha_;
  typename XViewType::const_type x_;
  typename YViewType::const_type y_;
  AViewType A_;
};

// Team parallel version of SYR2.
template <class ExecutionSpace, class XViewType, class YViewType, class AViewType, class IndexType, bool tJustTranspose,
          bool tJustUp>
void teamParallelSyr2(const ExecutionSpace& space, const typename AViewType::const_value_type& alpha,
                      const XViewType& x, const YViewType& y, const AViewType& A) {
  static_assert(std::is_integral<IndexType>::value, "IndexType must be an integer");

  using AlphaCoeffType = typename AViewType::non_const_value_type;

  if (x.extent(0) == 0) {
    // no entries to update
    return;
  } else if (y.extent(0) == 0) {
    // no entries to update
    return;
  } else if (alpha == Kokkos::ArithTraits<AlphaCoeffType>::zero()) {
    // no entries to update
    return;
  }

  constexpr bool isLayoutLeft = std::is_same<typename AViewType::array_layout, Kokkos::LayoutLeft>::value;
  using layout_tag =
      typename std::conditional<isLayoutLeft, TeamParallelSYR2_LayoutLeftTag, TeamParallelSYR2_LayoutRightTag>::type;
  using TeamPolicyType = Kokkos::TeamPolicy<ExecutionSpace, layout_tag>;
  TeamPolicyType teamPolicy;
  if (isLayoutLeft) {
    // LayoutLeft: one team per column
    teamPolicy = TeamPolicyType(space, A.extent(1), Kokkos::AUTO);
  } else {
    // LayoutRight: one team per row
    teamPolicy = TeamPolicyType(space, A.extent(0), Kokkos::AUTO);
  }

  TeamParallelSYR2<ExecutionSpace, XViewType, YViewType, AViewType, IndexType, tJustTranspose, tJustUp> functor(
      alpha, x, y, A);
  Kokkos::parallel_for("KokkosBlas::syr2[teamParallel]", teamPolicy, functor);
}

// ---------------------------------------------------------------------------------------------

// generalSyr2Impl():
// - use thread parallel code (rangePolicy) if execution space is CPU;
// - use team parallel code (teamPolicy) if execution space is GPU.
//
// The 'enable_if' makes sure unused kernels are not instantiated.

template <class ExecutionSpace, class XViewType, class YViewType, class AViewType, class IndexType, bool tJustTranspose,
          bool tJustUp,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>()>::type* = nullptr>
void generalSyr2Impl(const ExecutionSpace& space, const typename AViewType::const_value_type& alpha, const XViewType& x,
                     const YViewType& y, const AViewType& A) {
  threadParallelSyr2<ExecutionSpace, XViewType, YViewType, AViewType, IndexType, tJustTranspose, tJustUp>(space, alpha,
                                                                                                          x, y, A);
}

template <class ExecutionSpace, class XViewType, class YViewType, class AViewType, class IndexType, bool tJustTranspose,
          bool tJustUp,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>()>::type* = nullptr>
void generalSyr2Impl(const ExecutionSpace& space, const typename AViewType::const_value_type& alpha, const XViewType& x,
                     const YViewType& y, const AViewType& A) {
  teamParallelSyr2<ExecutionSpace, XViewType, YViewType, AViewType, IndexType, tJustTranspose, tJustUp>(space, alpha, x,
                                                                                                        y, A);
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS2_SYR2_IMPL_HPP_
