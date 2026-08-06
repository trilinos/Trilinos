// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef _KOKKOSGRAPH_Merge_IMPL_HPP
#define _KOKKOSGRAPH_Merge_IMPL_HPP

/*! \file KokkosGraph_Merge_impl.hpp
 *
 */

#include "KokkosGraph_MergePath.hpp"

namespace KokkosGraph {
namespace Impl {

/*! \brief sequential merge of sorted `a` and `b` into `c`

   No more than c.size() entries will be taken from `a` and `b` combined.
   If `a` or `b` are exhausted, fewer may be taken
   The number of entries taken from a or b depends on the values of those
   entries
*/
template <typename CView, typename AView, typename BView>
KOKKOS_INLINE_FUNCTION void merge_into_thread(const CView &c, const AView &a, const BView &b) {
  static_assert(AView::rank == 1,
                "merge_into_thread requires rank-1 A");  // .size()
  static_assert(BView::rank == 1, "merge_into_thread requires rank-1 B");
  static_assert(CView::rank == 1, "merge_into_thread requires rank-1 C");

  typename AView::size_type ai = 0;
  typename BView::size_type bi = 0;
  typename CView::size_type ci = 0;
  while (ai < a.size() && bi < b.size() && ci < c.size()) {
    const auto av = a(ai);
    const auto bv = b(bi);

    if (av < bv) {
      c(ci++) = av;
      ++ai;
    } else {
      c(ci++) = bv;
      ++bi;
    }
  }
  // a or b has been exhaused. merge from the other
  while (ai < a.size() && ci < c.size()) {
    c(ci++) = a(ai++);
  }
  while (bi < b.size() && ci < c.size()) {
    c(ci++) = b(bi++);
  }
}

template <typename TeamMember, typename CView, typename AView, typename BView>
KOKKOS_INLINE_FUNCTION void merge_into_stepper_team(const TeamMember &handle, const CView &c, const AView &a,
                                                    const BView &b, const size_t threadPathLength) {
  static_assert(AView::rank == 1,
                "merge_into_stepper_team requires rank-1 A");  // .size()
  static_assert(BView::rank == 1, "merge_into_team requires rank-1 B");
  static_assert(CView::rank == 1, "merge_into_team requires rank-1 C");

  using c_size_type = typename CView::size_type;
  using a_size_type = typename AView::size_type;
  using b_size_type = typename BView::size_type;

  auto stepper = [&](const StepDirection &dir, const StepperContext &step, const StepperContext &thread) {
    if (dir == StepDirection::a) {
      a_size_type ai = thread.ai + step.ai;
      c_size_type ci = thread.pi + step.pi;
      c(ci)          = a(ai);
    } else if (dir == StepDirection::b) {
      b_size_type bi = thread.bi + step.bi;
      c_size_type ci = thread.pi + step.pi;
      c(ci)          = b(bi);
    }
  };

  KokkosGraph::merge_path_team(handle, a, b, c.size(), stepper, threadPathLength);
}

template <typename TeamMember, typename CView, typename AView, typename BView>
KOKKOS_INLINE_FUNCTION void merge_into_raw_team(const TeamMember &handle, const CView &c, const AView &a,
                                                const BView &b, const size_t threadPathLength) {
  using c_size_type = typename CView::size_type;
  using a_size_type = typename AView::size_type;
  using b_size_type = typename BView::size_type;

  static_assert(AView::rank == 1,
                "merge_into_raw_team requires rank-1 A");  // .size()
  static_assert(BView::rank == 1, "merge_into_team requires rank-1 B");
  static_assert(CView::rank == 1, "merge_into_team requires rank-1 C");

  Kokkos::parallel_for(Kokkos::TeamThreadRange(handle, 0, handle.team_size()), [&](const size_t i) {
    const c_size_type diagonal = i * threadPathLength;

    if (diagonal >= c.size()) {
      return;  // not participating
    }

    // do diagonal search for work division
    auto dsr = Impl::diagonal_search(a, b, diagonal);

    // determine how long the thread arrays really are
    const c_size_type nC = KOKKOSKERNELS_MACRO_MIN(threadPathLength, c.size() - diagonal);
    const a_size_type nA = KOKKOSKERNELS_MACRO_MIN(threadPathLength, a.size() - dsr.ai);
    const b_size_type nB = KOKKOSKERNELS_MACRO_MIN(threadPathLength, b.size() - dsr.bi);

    auto cs = Kokkos::subview(c, Kokkos::make_pair(diagonal, diagonal + nC));
    auto as = Kokkos::subview(a, Kokkos::make_pair(dsr.ai, dsr.ai + nA));
    auto bs = Kokkos::subview(b, Kokkos::make_pair(dsr.bi, dsr.bi + nB));

    // sequential merge
    merge_into_thread(cs, as, bs);
  });
}

template <typename TeamMember, typename CView, typename AView, typename BView>
KOKKOS_INLINE_FUNCTION void merge_into_team(const TeamMember &handle, const CView &c, const AView &a, const BView &b,
                                            const size_t threadPathLength) {
  merge_into_raw_team(handle, c, a, b, threadPathLength);
}

/* merge a and b into c
 */
template <typename CView, typename AView, typename BView>
struct FlatMerger {
  constexpr FlatMerger(size_t chunkSize, const CView &c, const AView &a, const BView &b)
      : chunkSize_(chunkSize), c_(c), a_(a), b_(b) {}

  KOKKOS_INLINE_FUNCTION void operator()(const size_t i /*ith diagonal*/) const {
    const size_t diagonal = i * chunkSize_;
    if (diagonal >= c_.size()) {
      return;
    }

    // merge chunkSize values from the discovered starting place in a and b
    // where to start comparison in a and b
    auto dsr  = Impl::diagonal_search(a_, b_, diagonal);
    size_t ai = dsr.ai;
    size_t bi = dsr.bi;

    // unmanaged views for sequential merge
    // may not necessarily merge a whole chunk at the end of the merge list
    size_t nC = KOKKOSKERNELS_MACRO_MIN(chunkSize_, c_.size() - diagonal);
    // may may not be able to read from a or b due to their sizes
    size_t nA = KOKKOSKERNELS_MACRO_MIN(chunkSize_, a_.size() - ai);
    size_t nB = KOKKOSKERNELS_MACRO_MIN(chunkSize_, b_.size() - bi);

    auto sc = Kokkos::subview(c_, Kokkos::pair{diagonal, diagonal + nC});
    auto sa = Kokkos::subview(a_, Kokkos::pair{ai, ai + nA});
    auto sb = Kokkos::subview(b_, Kokkos::pair{bi, bi + nB});
    merge_into_thread(sc, sa, sb);
  }

  template <typename ExecSpace>
  static void launch(const ExecSpace &space, const CView &c, const AView &a, const BView &b) {
    // one diagonal for each unit of concurrency
    const size_t chunkSize = (c.size() + space.concurrency() - 1) / space.concurrency();
    const size_t numChunks = (c.size() + chunkSize - 1) / chunkSize;

    // Do the parallel merge
    static_assert(Kokkos::SpaceAccessibility<ExecSpace, typename AView::memory_space>::accessible,
                  "KokkoGraph::merge_into ExecSpace must be able to access A");
    static_assert(Kokkos::SpaceAccessibility<ExecSpace, typename BView::memory_space>::accessible,
                  "KokkoGraph::merge_into ExecSpace must be able to access B");
    static_assert(Kokkos::SpaceAccessibility<ExecSpace, typename CView::memory_space>::accessible,
                  "KokkoGraph::merge_into ExecSpace must be able to access C");
    Kokkos::RangePolicy<ExecSpace> policy(space, 0, numChunks);
    Impl::FlatMerger<CView, AView, BView> merger(chunkSize, c, a, b);
    Kokkos::parallel_for("KokkosGraph::merge_into", policy, merger);
  }

  size_t chunkSize_;
  CView c_;
  AView a_;
  BView b_;
};

template <typename ExecSpace, typename CView, typename AView, typename BView>
struct HierarchicalMerger {
  constexpr HierarchicalMerger(const CView &c, const AView &a, const BView &b, const size_t teamPathLength,
                               const size_t threadPathLength)
      : c_(c), a_(a), b_(b), teamPathLength_(teamPathLength), threadPathLength_(threadPathLength) {}

  using merger_type = HierarchicalMerger<ExecSpace, CView, AView, BView>;

  using execution_space    = ExecSpace;
  using scratch_space_type = typename execution_space::scratch_memory_space;

  using c_value_type = typename CView::non_const_value_type;
  using a_value_type = typename AView::non_const_value_type;
  using b_value_type = typename BView::non_const_value_type;

  using team_policy_type = Kokkos::TeamPolicy<execution_space>;
  using team_member_type = typename team_policy_type::member_type;

  using c_scratch_view = Kokkos::View<typename CView::data_type, scratch_space_type>;
  using a_scratch_view = Kokkos::View<typename AView::data_type, scratch_space_type>;
  using b_scratch_view = Kokkos::View<typename BView::data_type, scratch_space_type>;

  KOKKOS_INLINE_FUNCTION void operator()(const team_member_type &handle) const {
    // determine team's position in the merge path
    size_t teamDiagonal = handle.league_rank() * teamPathLength_;
    if (teamDiagonal >= c_.size()) {
      return;
    }
    auto dsr = Impl::diagonal_search(a_, b_, teamDiagonal);

    // figure out how much of view actually remain
    size_t nC = KOKKOSKERNELS_MACRO_MIN(teamPathLength_, c_.size() - teamDiagonal);
    size_t nA = KOKKOSKERNELS_MACRO_MIN(teamPathLength_, a_.size() - dsr.ai);
    size_t nB = KOKKOSKERNELS_MACRO_MIN(teamPathLength_, b_.size() - dsr.bi);

    // The entire path length may not be available since the
    // view dimensions may not divide evenly into teamPathLength_
    c_scratch_view cs(handle.team_scratch(0), nC);
    a_scratch_view as(handle.team_scratch(0), nA);
    b_scratch_view bs(handle.team_scratch(0), nB);

    // fill scratch
    Kokkos::parallel_for(Kokkos::TeamThreadRange(handle, 0, KOKKOSKERNELS_MACRO_MAX(as.size(), bs.size())),
                         [&](const size_t i) {
                           if (i < as.size()) {
                             as(i) = a_(dsr.ai + i);
                           }
                           if (i < bs.size()) {
                             bs(i) = b_(dsr.bi + i);
                           }
                         });
    handle.team_barrier();

    merge_into_team(handle, cs, as, bs, threadPathLength_);
    handle.team_barrier();

    // write back result
    Kokkos::parallel_for(Kokkos::TeamThreadRange(handle, 0, cs.size()),
                         [&](const size_t i) { c_(teamDiagonal + i) = cs(i); });
  }

  size_t team_shmem_size(int /*teamSize*/) const {
    return (sizeof(c_value_type) + sizeof(a_value_type) * sizeof(b_value_type)) * teamPathLength_;
  }

  static void launch(const ExecSpace &space, const CView &c, const AView &a, const BView &b) {
    static_assert(Kokkos::SpaceAccessibility<ExecSpace, typename AView::memory_space>::accessible,
                  "KokkoGraph::merge_into ExecSpace must be able to access A");
    static_assert(Kokkos::SpaceAccessibility<ExecSpace, typename BView::memory_space>::accessible,
                  "KokkoGraph::merge_into ExecSpace must be able to access B");
    static_assert(Kokkos::SpaceAccessibility<ExecSpace, typename CView::memory_space>::accessible,
                  "KokkoGraph::merge_into ExecSpace must be able to access C");

    // choose the team path length based on resources
    const size_t teamPathLength   = 256;
    const size_t leagueSize       = (c.size() + teamPathLength - 1) / teamPathLength;
    const size_t teamSize         = 128;
    const size_t threadPathLength = (teamPathLength + teamSize - 1) / teamSize;

    // Do the parallel merge
    team_policy_type policy(space, leagueSize, teamSize);
    merger_type merger(c, a, b, teamPathLength, threadPathLength);
    Kokkos::parallel_for("KokkosGraph::merge_into", policy, merger);
  }

  CView c_;
  AView a_;
  BView b_;
  size_t teamPathLength_;
  size_t threadPathLength_;
};

}  // namespace Impl
}  // namespace KokkosGraph

#endif  // _KOKKOSGRAPH_Merge_IMPL_HPP
