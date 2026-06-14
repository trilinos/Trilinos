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

#ifndef _KOKKOSGRAPH_MERGEPATH_HPP
#define _KOKKOSGRAPH_MERGEPATH_HPP

#include "KokkosKernels_SimpleUtils.hpp"
#include "KokkosGraph_MergePath_impl.hpp"

/*! \file KokkosGraph_MergePath.hpp
    \brief Provides a MergePath abstraction and related operations

    The "merge matrix" M of two sorted lists A and B has M[i,j] = 1 iff A[i] >
   B[j], and 0 otherwise.

    \verbatim
       0 1 2 3 B
      ________
    2| 1 1 0 0
    2| 1 1 0 0
    2| 1 1 0 0
    2| 1 1 0 0
    4| 1 1 1 1
    A
    \endverbatim

    The merge path follows the boundaries between the 0s and 1s.

    \verbatim
       0 1 2 3 B
      ________
      *->->
    2| 1 1|0 0
          v
    2| 1 1|0 0
          v
    2| 1 1|0 0
          v
    2| 1 1|0 0
          v->->
    4| 1 1 1 1|
    A         v
    \endverbatim

    This is equivalent to answering the question, "if I have two sorted lists A
   and B, what order should I take elements from them to produce a merged sorted
   list C", where ties are arbitrarily broken by choosing from list A.

    The length of the merge path is len(A) + len(B)

    This file provides functions that abstract over the merge path, calling a \c
   Stepper at each step along the merge path. The stepper is provided with a
   direction, as well as an index into list A, list B, and the number of steps
   along the merge path

    For the merge path shown above, the following Stepper calls are made:

    \verbatim
    stepper(b, {0,0,0})
    stepper(b, {0,1,1}) // step in b direction, +1 b, +1 path index
    stepper(a, {1,1,2}) // step in a direction, +1 a, +1 path index
    stepper(a, {2,1,3})
    stepper(a, {3,1,4})
    stepper(a, {4,1,5})
    stepper(b, {4,2,6})
    stepper(b, {4,3,7})
    stepper(a, {5,3,8}) // path length is 9, 9 steps have been made
    \endverbatim

    In practice, the indices into A, B, and the merge path are provided through
   \c StepperContext structs, and the direction through a \c StepDirection enum.

    The merge path is traversed hierarchically: each Team traverses a chunk, and
   each thread within the team a small part of the team's chunk. Therefore, the
   \c Stepper takes two contexts, a `step` and a `thread` context. The `step`
   context tells you where within the thread's piece the current step is. The
   `thread` context tells you where within the team's piece the current thread
   is. The global location in A, B, and the path can be recovered by summing the
   corresponding \c StepperContext members.
*/

namespace KokkosGraph {

/*! \brief Provides context for where in the merge matrix a step is taking place
 */
struct StepperContext {
  using offset_type = size_t;

  offset_type ai;  //!< index into a
  offset_type bi;  //!< index into b
  offset_type pi;  //!< index into the path

  KOKKOS_INLINE_FUNCTION
  constexpr StepperContext(offset_type _ai, offset_type _bi, offset_type _pi) : ai(_ai), bi(_bi), pi(_pi) {}

  KOKKOS_INLINE_FUNCTION
  constexpr StepperContext() : StepperContext(0, 0, 0) {}
};

/*! \enum StepDirection
    \brief Whether the merge-path step was along the A or B list
*/
enum class StepDirection { a, b };

/*! \brief Follow a merge path for two sorted input views and apply a stepper
   function to each element. \tparam AView Type of the first input view. \tparam
   BView Type of the second input view. \tparam Stepper Type of the stepper
   function. \tparam Ctxs Types of the optional contexts passed to the stepper
   function. \param a First sorted input view. \param b Second sorted input
   view. \param pathLength Maximum merge path length to process \param stepper
   Function to call for each path segment \param ctxs Optional contexts to pass
   to the stepper function.

    Follow a merge path for two sorted input views a and b and applies a stepper
   function at most pathLength elements. The stepper function should be
   invokable with a StepDirection enum value, followed by 0, 1, 2 StepperContext
   arguments.

    This generates the `step` StepperContext for each step of the merge path.
    The stepper will be invoked as stepper(StepDirection, step, ctx...)
*/
template <typename AView, typename BView, typename Stepper, typename... Ctxs>
KOKKOS_INLINE_FUNCTION void merge_path_thread(const AView &a, const BView &b, const size_t pathLength,
                                              const Stepper &stepper, Ctxs... ctxs) {
  static_assert(AView::rank == 1, "follow_path requires rank-1 AView");
  static_assert(BView::rank == 1, "follow_path requires rank-1 BView");

  constexpr bool ONE_CTX = std::is_invocable<Stepper, StepDirection, StepperContext>::value;
  constexpr bool TWO_CTX = std::is_invocable<Stepper, StepDirection, StepperContext, StepperContext>::value;
  constexpr bool THREE_CTX =
      std::is_invocable<Stepper, StepDirection, StepperContext, StepperContext, StepperContext>::value;
  static_assert(ONE_CTX || TWO_CTX || THREE_CTX,
                "Stepper should be invokable with a StepDirection, then 1, 2, "
                "or 3 StepperContext arguments, for the step, thread, and team "
                "context respectively");

  static_assert(sizeof...(Ctxs) == 0 || sizeof...(Ctxs) == 1 || sizeof...(Ctxs) == 2,
                "Zero, one, or two contexts should be passed to merge_path_thread");

  StepperContext step;
  while (step.pi < pathLength && step.ai < a.size() && step.bi < b.size()) {
    if (a(step.ai) <= b(step.bi)) {  // step in A direction
      stepper(StepDirection::a, step, ctxs...);
      ++step.ai;
      ++step.pi;
    } else {  // step in B direction
      stepper(StepDirection::b, step, ctxs...);
      ++step.bi;
      ++step.pi;
    }
  }
  while (step.pi < pathLength && step.ai < a.size()) {
    stepper(StepDirection::a, step, ctxs...);
    ++step.ai;
    ++step.pi;
  }
  while (step.pi < pathLength && step.bi < b.size()) {
    stepper(StepDirection::b, step, ctxs...);
    ++step.bi;
    ++step.pi;
  }
}

/*! \brief Collaboratively follow a merge path for two sorted input views and
   apply a stepper function to each element. \tparam AView A Kokkos::View
    \tparam BViewLike A Kokkos::View or KokkosKernels::Impl::Iota
    \tparam Stepper Type of the stepper function.
    \tparam Ctxs Types of the optional contexts passed to the stepper function.
    \param a First sorted input view.
    \param b Second sorted input view.
    \param pathLength Maximum merge path length to process
    \param stepper Function to call for each path segment
    \param ctxs Optional contexts to pass to the stepper function.
    \param threadPathLength Maximum length each thread will process

    Collaboartively calls merge_path_thread on subsegments, adding a `thread`
   context to the stepper call that says where along the merge path that
   thread's processing begins
*/
template <typename TeamHandle, typename AView, typename BViewLike, typename Stepper>
KOKKOS_INLINE_FUNCTION void merge_path_team(const TeamHandle &handle, const AView &a, const BViewLike &b,
                                            const size_t pathLength, const Stepper &stepper, size_t threadPathLength) {
  static_assert(AView::rank == 1, "merge_path_team requires rank-1 a");
  static_assert(BViewLike::rank == 1, "merge_path_team requires rank-1 b");

  Kokkos::parallel_for(Kokkos::TeamThreadRange(handle, 0, handle.team_size()), [&](const size_t i) {
    // split up with a diagonal search
    // size_t threadPathLength =
    //     (pathLength + handle.team_size() - 1) / handle.team_size();
    const size_t diagonal = threadPathLength * i;
    if (diagonal >= pathLength) {
      return;
    }
    auto dsr  = Impl::diagonal_search(a, b, diagonal);
    size_t ai = dsr.ai;
    size_t bi = dsr.bi;

    // capture where in the team context the thread is working
    const StepperContext threadCtx(ai, bi, diagonal);

    // final piece pay be shorter
    threadPathLength = KOKKOSKERNELS_MACRO_MIN(threadPathLength, pathLength - diagonal);

    // take appropriate subviews of A and B according to diagonal search
    size_t nA = KOKKOSKERNELS_MACRO_MIN(threadPathLength, a.size() - ai);
    size_t nB = KOKKOSKERNELS_MACRO_MIN(threadPathLength, b.size() - bi);
    auto as   = Kokkos::subview(a, Kokkos::make_pair(ai, ai + nA));
    if constexpr (KokkosKernels::Impl::is_iota<BViewLike>::value) {
      BViewLike bs(b, Kokkos::make_pair(bi, bi + nB));
      merge_path_thread(as, bs, threadPathLength, stepper, threadCtx);
    } else {
      auto bs = Kokkos::subview(b, Kokkos::make_pair(bi, bi + nB));

      // each thread contributes a path segment in parallel
      merge_path_thread(as, bs, threadPathLength, stepper, threadCtx);
    }
  });
}

template <typename TeamHandle, typename AView, typename BViewLike, typename Stepper>
KOKKOS_INLINE_FUNCTION void merge_path_team(const TeamHandle &handle, const AView &a, const BViewLike &b,
                                            const size_t pathLength, const Stepper &stepper) {
  const size_t threadPathLength = (pathLength + handle.team_size() - 1) / handle.team_size();
  merge_path_team(handle, a, b, pathLength, stepper, threadPathLength);
}

}  // namespace KokkosGraph

#endif  // _KOKKOSGRAPH_MERGEPATH_HPP
