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

#ifndef _KOKKOSGRAPH_LOADBALANCE_HPP
#define _KOKKOSGRAPH_LOADBALANCE_HPP

#include <Kokkos_Core.hpp>

#include "KokkosKernels_SimpleUtils.hpp"

#include "KokkosGraph_LoadBalance_impl.hpp"
#include "KokkosGraph_MergePath.hpp"

/*! \file KokkosGraph_LoadBalance.hpp
    \brief Implementation of global, team, and thread load-balance

    Consider a scenario where you have a list of irregularly-sized tasks; the
   0th task has size 2 (2 work items), the first task has size 1, the second 0,
   and the third 3.

    \verbatim
    taskSizes = {2, 1, 0, 3}
    \endverbatim
    This is a total of 6 = 2 + 1 + 0 + 3 work-items.

    Let's say for each work-item you want to know which task it came from, and
   what work-item number it is within that task. For example, the first two
   work-items are the 0th and 1st work-items in task 0. The next work item is
   the 0th work item in task 1, and the final three work-items are numbers 0, 1,
   and 2 in task 3.

    \verbatim
    tasks = {0,0,1,3,3,3}
    ranks = {0,1,0,0,1,2}
    \endverbatim

    This file provides global, team, and thread interfaces for functions to do
   that.
*/

namespace KokkosGraph {

/*! \struct LoadBalanceResult A load-balance function result, where work-item
   `i` came from Task = \c tasks(i) and was the `ranks(i)`-th work-item in that
   task


   The result of calling load-balance on a TaskSizesView of n tasks, where each
   entry in the view corresponds to the number of work-items (or size) of the
   task, `tasks(i)` and `ranks(i)` are the source task for work-item, and which
   work item it was within the task, respectively.

   For example, if the task-sizes view was {2, 1, 0, 3}
   tasks = {0, 0, 1, 3, 3, 3} and ranks = {0, 1, 0, 0, 1, 2}
   i.e., work-item 4 has index 1 in task 3

   \tparam TaskSizesView The input task-sizes view, which determines the struct
   view types

   TODO: add labels to views
*/
template <typename TaskSizesView>
struct LoadBalanceResult {
  using tasks_view_type = typename Impl::ResultTypes<TaskSizesView>::tasks_view_type;
  using ranks_view_type = typename Impl::ResultTypes<TaskSizesView>::ranks_view_type;

  tasks_view_type tasks;  ///< the task each work item came from
  ranks_view_type ranks;  ///< the rank of the work-item within that task

  /*! \brief resize the members of this struct to support totalWork work-items
   */
  KOKKOS_INLINE_FUNCTION
  void resize_for_total_work(const typename TaskSizesView::non_const_value_type &totalWork) {
    Kokkos::resize(tasks, totalWork);
    Kokkos::resize(ranks, totalWork);
  }
};

/*! \brief Produce a LoadBalanceResult on an exclusive prefix sum of the task
   sizes

    \tparam ExecSpace Type of the execution space instance
    \tparam View Kokkos::View of task sizes

    \param[in] space Execution space instance to work in
    \param[in] scanTasks exclusive prefix sum of the task sizes
    \param[in] totalWork sum of the task sizes
 */
template <typename ExecSpace, typename View>
LoadBalanceResult<View> load_balance_exclusive(const ExecSpace &space, const View &scanTasks,
                                               const typename View::non_const_value_type totalWork) {
  static_assert(View::rank == 1, "scanTasks must be rank 1");

  // don't include the leading zero
  auto tail = Kokkos::subview(scanTasks, Kokkos::make_pair(size_t(1), scanTasks.size() - 1));
  LoadBalanceResult<View> lbr;
  Impl::load_balance_inclusive(space, lbr.tasks, lbr.ranks, tail, totalWork);
  return lbr;
}

/*! \brief Produce a LoadBalanceResult on an exclusive prefix sum of the task
   sizes

    \tparam View Kokkos::View of task sizes

    \param[in] scanTasks exclusive prefix sum of the task sizes
    \param[in] totalWork sum of the task sizes
 */
template <typename View>
LoadBalanceResult<View> load_balance_exclusive(const View &scanTasks,
                                               const typename View::non_const_value_type totalWork) {
  using execution_space = typename View::device_type::execution_space;
  return load_balance_exclusive(execution_space(), scanTasks, totalWork);
}

/*! \brief Produce a LoadBalanceResult from a view of task sizes

    \tparam ExecSpace Type of the execution space instance
    \tparam View Kokkos::View of task sizes

    \param[in] space Execution space instance to work in
    \param[in] taskSizes the number of work-items in each task
    \param[in] totalWork sum of the task sizes
 */
template <typename ExecSpace, typename View>
LoadBalanceResult<View> load_balance(const ExecSpace &space, const View &taskSizes) {
  using value_type = typename View::non_const_value_type;

  static_assert(View::rank == 1, "KokkoGraph::load_balance requires rank-1 taskSizes");
  static_assert(Kokkos::SpaceAccessibility<ExecSpace, typename View::memory_space>::accessible,
                "ExecSpace must be able to access taskSizes");

  if (0 == taskSizes.size()) {
    return LoadBalanceResult<View>();
  }

  View sum("task-sizes-prefix-sum", taskSizes.size());
  Kokkos::deep_copy(space, sum, taskSizes);
  KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum(space, sum.size(), sum);

  // retrieve work items sum from view's memory space
  value_type numWorkItems;
  Kokkos::deep_copy(space, numWorkItems, Kokkos::subview(sum, sum.size() - 1));
  space.fence();

  LoadBalanceResult<View> lbr;
  Impl::load_balance_inclusive(space, lbr.tasks, lbr.ranks, sum, numWorkItems);
  return lbr;
}

/*! \brief Produce a LoadBalanceResult from a view of task sizes

    \tparam View Kokkos::View of task sizes

    \param[in] taskSizes the number of work-items in each task
 */
template <typename View>
LoadBalanceResult<View> load_balance(const View &taskSizes) {
  typedef typename View::device_type::execution_space execution_space;
  return load_balance(execution_space(), taskSizes);
}

template <typename View>
LoadBalanceResult<View> load_balance_thread(const View &taskSizes) {
  // Each entry in the view represents a work item, thus the size_type of the
  // view is how the work items are numbered The value of the entry is the
  // amount of work in that work item
  using work_item_size = typename View::non_const_value_type;

  // to be returned
  LoadBalanceResult<View> lbr;

  // get total work size
  work_item_size totalWork = 0;
  for (size_t i = 0; i < taskSizes.size(); ++i) {
    totalWork += taskSizes(i);
  }

  // size result appropriately
  lbr.resize_for_total_work(totalWork);

  // do the load-balancing
  Impl::load_balance_into_thread(lbr.assignment, lbr.rank, taskSizes);
  return lbr;
}

/*! \brief team-collaborative load-balance into preallocated views

    \tparam TeamMember Type of `handle`
    \tparam TaskView type of `tasks` view
    \tparam RanksView type of `ranks` view
    \tparam TaskSizesView type of `inclusivePrefixSumTaskSizes`

    \param[in] handle the Kokkos::TeamPolicy handle
    \param[out] tasks entry i contains the source task for work-item i
    \param[out] ranks entry i contains the rank of work-item i in it's source
   task \param[in] inclusivePrefixSumTaskSizes inclusive prefix sum of task
   sizes

    inclusive_prefix_sum_team may be used for a team-collaborative inclusive
   prefix sum to provide `inclusivePrefixSumTaskSizes`
*/
template <typename TeamMember, typename TasksView, typename RanksView, typename TaskSizesView>
KOKKOS_INLINE_FUNCTION void load_balance_into_team(const TeamMember &handle, const TasksView &tasks,
                                                   const RanksView &ranks,
                                                   const TaskSizesView &inclusivePrefixSumTaskSizes) {
  using task_size_type = typename TaskSizesView::non_const_value_type;
  using iota_type      = KokkosKernels::Impl::Iota<task_size_type>;

  if (inclusivePrefixSumTaskSizes.size() == 0) {
    return;
  }

  auto stepper = [&](const StepDirection &dir, const StepperContext &step, const StepperContext &thread) {
    if (StepDirection::b == dir) {
      // recover global index into a and b
      auto ag   = thread.ai + step.ai;
      auto bg   = thread.bi + step.bi;
      tasks(bg) = ag;
      ranks(bg) = bg - (ag > 0 ? inclusivePrefixSumTaskSizes(ag - 1) : 0);
    }
  };

  const task_size_type totalWork = inclusivePrefixSumTaskSizes(inclusivePrefixSumTaskSizes.size() - 1);
  merge_path_team(handle, inclusivePrefixSumTaskSizes, iota_type(totalWork),
                  inclusivePrefixSumTaskSizes.size() + totalWork, stepper);
}

/*! \brief team-collaborative inclusive prefix sum

    \tparam TeamMember Type of `handle`
    \tparam OutView type of `tasks` view
    \tparam InView type of `ranks` view

    \param[in] handle the Kokkos::TeamPolicy handle
    \param[out] out the inclusive prefix sum
    \param[in] in the values to be summed
*/
template <typename TeamMember, typename OutView, typename InView>
KOKKOS_INLINE_FUNCTION void inclusive_prefix_sum_team(const TeamMember &handle, const OutView &out, const InView &in) {
  using index_type = typename OutView::size_type;
  using value_type = typename OutView::non_const_value_type;

  Kokkos::parallel_scan(
      Kokkos::TeamThreadRange(handle, 0, out.size()),
      KOKKOS_LAMBDA(const index_type i, value_type &partial, const bool final) {
        partial += in(i);
        if (final) {
          out(i) = partial;
        }
      });
}

/*! \brief team-collaborative load-balance
    \tparam TeamMember Type of `handle`
    \tparam View Type of `taskSizes`

    \param[in] handle the Kokkos::TeamPolicy handle
    \param[in] taskSizes The number of work-items in each task
*/
template <typename TeamMember, typename View>
LoadBalanceResult<View> load_balance_team(const TeamMember &handle, const View &taskSizes) {
  using task_size_type = typename View::non_const_value_type;

  LoadBalanceResult<View> lbr;
  if (0 == taskSizes.size()) {
    return lbr;
  }
  View scratch("scratch", taskSizes.size());
  inclusive_prefix_sum_team(handle, scratch, taskSizes);
  const task_size_type totalWork = scratch(scratch.size() - 1);

  Kokkos::single(
      Kokkos::PerTeam(handle), KOKKOS_LAMBDA() { lbr.resize_for_total_work(totalWork); });
  load_balance_into_team(handle, lbr.tasks, lbr.ranks, scratch);
  return lbr;
}

}  // namespace KokkosGraph

#endif  // _KOKKOSGRAPH_LOADBALANCE_HPP
