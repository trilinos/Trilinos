// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef _KOKKOSGRAPH_LOADBALANCE_IMPL_HPP
#define _KOKKOSGRAPH_LOADBALANCE_IMPL_HPP

/*! \file KokkosGraph_LoadBalance_impl.hpp
 *
 */

#include "KokkosKernels_Iota.hpp"
#include "KokkosGraph_MergePath.hpp"
#include "KokkosKernels_SimpleUtils.hpp"

namespace KokkosGraph {
namespace Impl {

template <typename TaskSizesView>
struct ResultTypes {
  // each worker is assigned to an index of TaskSizesView
  using tasks_view_type = Kokkos::View<typename TaskSizesView::size_type *, typename TaskSizesView::device_type>;
  // each worker is assigned a rank within a work item, counted the same way
  // work items are
  using ranks_view_type =
      Kokkos::View<typename TaskSizesView::non_const_value_type *, typename TaskSizesView::device_type>;
};

/*! \brief

    \tparam scanTasks inclusive prefix sum of task sizes, so the final entry is
   the total work

    TODO: construct Iota on the fly in the body
*/
template <typename TaskView,      // type of assignment view
          typename RankView,      // type of rank view
          typename ScanTasksView  // type of lists to merge
          >
struct Balancer {
  Balancer(size_t chunkSize,
           const TaskView &tasks,  // output
           const RankView &ranks,  // output
           const ScanTasksView &scanTasks)
      : chunkSize_(chunkSize), tasks_(tasks), ranks_(ranks), scanTasks_(scanTasks) {}

  typedef typename ScanTasksView::size_type size_type;
  typedef typename ScanTasksView::non_const_value_type value_type;
  using iota_type = KokkosKernels::Impl::Iota<value_type, typename ScanTasksView::size_type>;

  KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const {
    typename iota_type::size_type totalWork = scanTasks_(scanTasks_.size() - 1);
    iota_type iota(totalWork);

    const size_t diagonal = i * chunkSize_;
    if (diagonal >= scanTasks_.size() + iota.size()) {
      return;
    }

    // diagonal search to figure out what part of the input the thread will
    // load-balance
    auto dsr  = Impl::diagonal_search(scanTasks_, iota, diagonal);
    size_t si = dsr.ai;
    size_t ii = dsr.bi;

#if 0
        size_t nS = KOKKOSKERNELS_MACRO_MIN(chunkSize_, scanTasks_.size() - si);
        size_t nI = KOKKOSKERNELS_MACRO_MIN(chunkSize_, iota.size() - ii);
        auto ss = Kokkos::subview(scanTasks_, Kokkos::pair{si, si+nS});
        auto is = iota_type(nI, ii);
        StepperContext threadCtx(si, ii, diagonal);

        /* We will provide the inclusive scan of task sizes as A, and Iota as B
        Therefore, stepping in B means there's another work-item in the current task
        
        This work-item comes from task at the global index into A
        This work-item is rank bg - scanTasks(ag-1)
        */
        auto stepper = [&](const StepDirection &dir,
                            const StepperContext &step,
                            const StepperContext &thread,
                            const StepperContext &team) {
            if (StepDirection::b == dir) {
            // recover global index into a and b
            size_t ag = team.ai + thread.ai + step.ai;
            size_t bg = team.bi + thread.bi + step.bi;
            tasks_(bg) = ag;
            ranks_(bg) = bg - (ag > 0 ? scanTasks_(ag - 1) : 0);
            }
        };


        merge_path_thread(ss, is, chunkSize_, NoopStepper(), bstepper, threadCtx);
#else
    // printf("d=%lu merge from si=%lu ii=%lu\n", diagonal, si, ii);
    // don't make more than chunk-size steps along the merge-path
    for (size_t m = 0; m < chunkSize_ && si < scanTasks_.size() && ii < iota.size();) {
      value_type sv = scanTasks_(si);
      value_type iv = iota(ii);

      // printf("(1) d=%lu si=%lu ii=%lu sv=%d iv=%d\n", diagonal, si, ii,
      // int(sv), int(iv));

      if (iv < sv) {
        // printf("(2) d=%lu ii=%lu tasks_.size()=%lu\n", diagonal, ii,
        // tasks_.size());
        tasks_(ii) = si;
        // printf("(3) d=%lu ii=%lu \n", diagonal, ii);
        ranks_(ii) = ii - (si > 0 ? scanTasks_(si - 1) : 0);  // implicit prefixed 0 on inclusive scan
        // printf("(4) d=%lu ii=%lu \n", diagonal, ii);
        ++ii;
        ++m;
      } else {
        ++si;
      }
    }
    // printf("(5) d=%lu\n", diagonal);
#endif
  }

  size_t chunkSize_;
  TaskView tasks_;
  RankView ranks_;
  ScanTasksView scanTasks_;
};

/*! \brief

    \param
*/
template <typename TaskView, typename RankView, typename View>
KOKKOS_INLINE_FUNCTION void load_balance_into_thread(const TaskView &tasks, const RankView &ranks,
                                                     const View &taskSizes) {
  using task_size_type = typename RankView::non_const_value_type;

  size_t i = 0;
  for (size_t task = 0; task < taskSizes.size(); ++task) {
    auto taskSize = taskSizes(task);

    // generate one piece of work for each size of the work item
    // each piece of work has a source work item and a rank within that work
    // item
    for (task_size_type rank = 0; rank < task_size_type(taskSize); ++rank) {
      tasks(i) = task;
      ranks(i) = rank;
      ++i;
    }
  }
}

/* `scanItems` is the inclusive prefix sum of the `workItems` view described for
 * `KokkosGraph::load_balance`
 */
template <typename ExecSpace, typename View>
void load_balance_inclusive(const ExecSpace &space, typename ResultTypes<View>::tasks_view_type &tasks,
                            typename ResultTypes<View>::ranks_view_type &ranks, const View &scanTasks,
                            const typename View::non_const_value_type totalWork) {
  static_assert(View::rank == 1, "scanTasks must be rank 1");
  static_assert(Kokkos::SpaceAccessibility<ExecSpace, typename View::memory_space>::accessible, "");
  static_assert(
      Kokkos::SpaceAccessibility<ExecSpace, typename ResultTypes<View>::tasks_view_type::memory_space>::accessible, "");
  static_assert(
      Kokkos::SpaceAccessibility<ExecSpace, typename ResultTypes<View>::ranks_view_type::memory_space>::accessible, "");

  Kokkos::resize(Kokkos::WithoutInitializing, tasks, totalWork);
  Kokkos::resize(Kokkos::WithoutInitializing, ranks, totalWork);
  if (0 == totalWork) {
    return;
  }

  const size_t chunkSize = (totalWork + space.concurrency() - 1) / space.concurrency();
  Balancer balancer(chunkSize, tasks, ranks, scanTasks);
  Kokkos::RangePolicy<ExecSpace> policy(space, 0, (scanTasks.size() + totalWork + chunkSize - 1) / chunkSize);
  Kokkos::parallel_for("KokkosGraph::load_balance_inclusive", policy, balancer);
}

}  // namespace Impl
}  // namespace KokkosGraph

#endif  // _KOKKOSGRAPH_LOADBALANCE_IMPL_HPP
