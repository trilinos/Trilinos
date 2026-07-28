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

#ifndef _KOKKOSGRAPH_MERGE_HPP
#define _KOKKOSGRAPH_MERGE_HPP

#include <Kokkos_Core.hpp>

#include "KokkosKernels_SimpleUtils.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"

#include "KokkosGraph_Merge_impl.hpp"
#include "KokkosGraph_MergePath.hpp"

namespace KokkosGraph {

/*! \file KokkosGraph_Merge.hpp
    \brief Implementation of global, team, and thread merge of sorted lists

    In general, these functions all operate on two lists sorted by the less-than
   operation, i.e. a[i] !< a[j] for all i > j They produce a merged list sorted
   in the same way. The entries within the merged result are stable with respect
   to the input lists.

    \verbatim
    a   = {0   2 3 5}
    b   = {  1   3  }
    out = {0 1 2 3 3 5}
    \endverbatim
*/

/*! \brief Sequential merge of sorted `a` and `b` into `c`

    \tparam CView Rank-1 Kokkos::View
    \tparam AView Rank-1 Kokkos::View
    \tparam BView Rank-1 Kokkos::View
    \param[out] c resulting merged \c a and \c b
    \param[in] a sorted view to merge
    \param[in] b sorted view to merge

   No more than c.size() entries will be taken from `a` and `b` combined.
   If `a` or `b` are exhausted, fewer may be taken
   The number of entries taken from a or b depends on the values of those
   entries
*/
template <typename CView, typename AView, typename BView>
KOKKOS_INLINE_FUNCTION void merge_into_thread(const CView &c, const AView &a, const BView &b) {
  Impl::merge_into_thread(c, a, b);
}

/*! \brief Team-collaborative parallel merge of sorted `a` and `b` into `c`

    \tparam TeamMember the Kokkos::TeamPolicy member type
    \tparam CView Rank-1 Kokkos::View
    \tparam AView Rank-1 Kokkos::View
    \tparam BView Rank-1 Kokkos::View
    \param[in] handle the team policy member
    \param[out] c resulting merged \c a and \c b
    \param[in] a sorted view to merge
    \param[in] b sorted view to merge
    \param[in] threadPathLength ceil((a.size() + b.size()) / handle.team_size())

   No more than c.size() entries will be taken from `a` and `b` combined.
   If `a` or `b` are exhausted, fewer may be taken
   The number of entries taken from a or b depends on the values of those
   entries

   \c threadPathLength may be provided as an optimization if it is precomputed
   If it is not known, consider a different overload
*/
template <typename TeamMember, typename CView, typename AView, typename BView>
KOKKOS_INLINE_FUNCTION void merge_into_team(const TeamMember &handle, const CView &c, const AView &a, const BView &b,
                                            const size_t threadPathLength) {
  return Impl::merge_into_team(handle, c, a, b, threadPathLength);
}

/*! \brief Team-collaborative parallel merge of sorted `a` and `b` into `c`

    \tparam TeamMember the Kokkos::TeamPolicy member type
    \tparam CView Rank-1 Kokkos::View
    \tparam AView Rank-1 Kokkos::View
    \tparam BView Rank-1 Kokkos::View
    \param[in] handle the team policy member
    \param[out] c resulting merged \c a and \c b
    \param[in] a sorted view to merge
    \param[in] b sorted view to merge

   No more than c.size() entries will be taken from `a` and `b` combined.
   If `a` or `b` are exhausted, fewer may be taken
   The number of entries taken from a or b depends on the values of those
   entries
*/
template <typename TeamMember, typename CView, typename AView, typename BView>
KOKKOS_INLINE_FUNCTION void merge_into_team(const TeamMember &handle, const CView &c, const AView &a, const BView &b) {
  using size_type = typename CView::size_type;
  static_assert(CView::rank == 1, "merge_into_team requires rank-1 C");
  const size_type threadPathLength = (c.size() + handle.team_size() - 1) / handle.team_size();
  merge_into_team(handle, c, a, b, threadPathLength);
}

template <typename ExecSpace, typename CView, typename AView, typename BView>
void merge_into(const ExecSpace &space, const CView &c, const AView &a, const BView &b) {
  static_assert(AView::rank == 1,
                "KokkosGraph::merge_into requires rank-1 A views");  // .size()
  static_assert(BView::rank == 1,
                "KokkosGraph::merge_into requires rank-1 B views");  // .size()
  static_assert(CView::rank == 1,
                "KokkosGraph::merge_into requires rank-1 C views");  // .size()

  if (c.size() != a.size() + b.size()) {
    throw std::logic_error("KokkosGraph::merge_into: c.size() must be a.size() + b.size()");
  }

  // if one or boths Views are empty, no real merging to do
  if (0 == a.size() && 0 == b.size()) {
    return;
  } else if (0 == a.size() && 0 != b.size()) {
    Kokkos::deep_copy(space, c, b);
    return;
  } else if (0 != a.size() && 0 == b.size()) {
    Kokkos::deep_copy(space, c, a);
    return;
  }

  // create unmanaged views
  using u_aview_t =
      Kokkos::View<typename AView::data_type, typename AView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using u_bview_t =
      Kokkos::View<typename BView::data_type, typename BView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using u_cview_t =
      Kokkos::View<typename CView::data_type, typename CView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  u_cview_t uc(c);
  u_aview_t ua(a);
  u_bview_t ub(b);

  if (KokkosKernels::Impl::is_gpu_exec_space_v<ExecSpace>) {
    Impl::HierarchicalMerger<ExecSpace, u_cview_t, u_aview_t, u_bview_t>::launch(space, uc, ua, ub);
  } else {
    Impl::FlatMerger<u_cview_t, u_aview_t, u_bview_t>::launch(space, uc, ua, ub);
  }
}

/*! \brief Global parallel merge of sorted `a` and `b` into `c`

    \tparam CView Rank-1 Kokkos::View
    \tparam AView Rank-1 Kokkos::View
    \tparam BView Rank-1 Kokkos::View
    \param[out] c resulting merged \c a and \c b
    \param[in] a sorted view to merge
    \param[in] b sorted view to merge

   No more than c.size() entries will be taken from `a` and `b` combined.
   If `a` or `b` are exhausted, fewer may be taken
   The number of entries taken from a or b depends on the values of those
   entries
*/
template <typename CView, typename AView, typename BView>
void merge_into(const CView &c, const AView &a, const BView &b) {
  using execution_space = typename AView::execution_space;
  merge_into(execution_space(), c, a, b);
}

/*! \brief Global merge sorted a and b, resizing c appropriately

    \tparam ExecSpace A Kokkos execution space
    \tparam CView Rank-1 Kokkos::View
    \tparam AView Rank-1 Kokkos::View
    \tparam BView Rank-1 Kokkos::View
    \param[in] space the execution space instance to work in
    \param[out] c resulting merged \c a and \c b
    \param[in] a sorted view to merge
    \param[in] b sorted view to merge

    \c c will be resized to a.size() + b.size() before the merge is invoked
*/
template <typename ExecSpace, typename CView, typename AView, typename BView>
void resize_and_merge_into(const ExecSpace &space, CView &c, const AView &a, const BView &b) {
  static_assert(AView::rank == 1, "KokkoGraph::merge requires rank-1 A");
  static_assert(BView::rank == 1, "KokkoGraph::merge requires rank-1 B");
  static_assert(CView::rank == 1, "KokkoGraph::merge requires rank-1 C");
  Kokkos::resize(c, a.size() + b.size());
  merge_into(space, c, a, b);
}

/*! \brief Global merge sorted a and b, resizing c appropriately

    \tparam CView Rank-1 Kokkos::View
    \tparam AView Rank-1 Kokkos::View
    \tparam BView Rank-1 Kokkos::View
    \param[out] c resulting merged \c a and \c b
    \param[in] a sorted view to merge
    \param[in] b sorted view to merge

    \c c will be resized to a.size() + b.size() before the merge is invoked
*/
template <typename CView, typename AView, typename BView>
void resize_and_merge_into(CView &c, const AView &a, const BView &b) {
  using execution_space = typename CView::device_type::execution_space;
  using u_aview_t =
      Kokkos::View<typename AView::data_type, typename AView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using u_bview_t =
      Kokkos::View<typename BView::data_type, typename BView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  u_aview_t ua(a);
  u_bview_t ub(b);
  resize_and_merge_into(execution_space(), c, ua, ub);
}

/*! \brief Global merge sorted a and b

    \tparam CView Rank-1 Kokkos::View
    \tparam ExecSpace Execution space used for the merge operation.
    \tparam AView Rank-1 Kokkos::View
    \tparam BView Rank-1 Kokkos::View
    \param[in] space the execution space instance to work in
    \param[in] a sorted view to merge
    \param[in] b sorted view to merge
    \return A CView of merged \c a and \c b, size a.size() + b.size()

    Merges two sorted input views a and b into a single sorted output view c.
    The output view is allocated and returned to the caller.
    The input views must be rank-1.
*/
template <typename CView, typename ExecSpace, typename AView, typename BView>
CView merge(const ExecSpace &space, const AView &a, const BView &b) {
  using u_aview_t =
      Kokkos::View<typename AView::data_type, typename AView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using u_bview_t =
      Kokkos::View<typename BView::data_type, typename BView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  u_aview_t ua(a);
  u_bview_t ub(b);

  static_assert(AView::rank == 1, "KokkoGraph::merge requires rank-1 A");
  CView c(Kokkos::ViewAllocateWithoutInitializing("merge_result"), a.size() + b.size());
  merge_into(space, c, ua, ub);
  return c;
}

/*! \brief Global merge sorted a and b

    \tparam CView Rank-1 Kokkos::View
    \tparam AView Rank-1 Kokkos::View
    \tparam BView Rank-1 Kokkos::View
    \param[in] a sorted view to merge
    \param[in] b sorted view to merge
    \return A CView of merged \c a and \c b, size a.size() + b.size()

    Merge two sorted input views a and b and returns a sorted view.
*/
template <typename CView, typename AView, typename BView>
CView merge(const AView &a, const BView &b) {
  using execution_space = typename AView::execution_space;
  using u_aview_t =
      Kokkos::View<typename AView::data_type, typename AView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using u_bview_t =
      Kokkos::View<typename BView::data_type, typename BView::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  u_aview_t ua(a);
  u_bview_t ub(b);
  return merge<CView>(execution_space(), ua, ub);
}

}  // namespace KokkosGraph

#endif  // _KOKKOSGRAPH_MERGE_HPP
