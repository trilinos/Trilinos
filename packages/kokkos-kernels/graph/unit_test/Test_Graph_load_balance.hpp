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

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include "KokkosGraph_LoadBalance.hpp"

#include <vector>

template <typename TaskSize, typename TaskIndex>
struct LBR {
  using task_index_type = TaskIndex;
  using task_size_type  = TaskSize;

  std::vector<task_index_type> tasks;
  std::vector<task_size_type> ranks;
};

template <typename TaskSize, typename TaskIndex>
LBR<TaskSize, TaskIndex> vector_load_balance(const std::vector<TaskSize> &taskSizes) {
  LBR<TaskSize, TaskIndex> lbr;

  for (TaskIndex task = 0; task < taskSizes.size(); ++task) {
    auto taskSize = taskSizes[task];

    // generate one piece of work for each size of the work item
    // each piece of work has a source work item and a rank within that work
    // item
    for (TaskSize rank = 0; rank < taskSize; ++rank) {
      lbr.tasks.push_back(task);
      lbr.ranks.push_back(rank);
    }
  }

  return lbr;
}

/* compare expected and actual
   if errors, print a, b, expected, and actual
*/
template <typename TaskSize, typename TaskIndex>
bool vector_compare(const std::vector<TaskSize> &taskSizes, const LBR<TaskSize, TaskIndex> &exp,
                    const LBR<TaskSize, TaskIndex> &act) {
  size_t err = 0;
  for (size_t i = 0; i < act.tasks.size() || i < act.ranks.size() || i < exp.tasks.size() || i < exp.ranks.size();
       ++i) {
    if (i < exp.tasks.size() && i < act.tasks.size()) {
      if (act.tasks[i] != exp.tasks[i]) {
        ++err;  // value mimatch
      }
    } else {
      ++err;  // size mismatch
    }
    if (i < exp.ranks.size() && i < act.ranks.size()) {
      if (act.ranks[i] != exp.ranks[i]) {
        ++err;  // value mimatch
      }
    } else {
      ++err;  // size mismatch
    }
  }
  if (err != 0) {
    std::cerr << "i\tts\te ta.\ta ta.\te ra.\ta ra." << std::endl;
    std::cerr << "-\t--\t-----\t-----\t-----\t-----" << std::endl;
    for (size_t i = 0; i < act.tasks.size() || i < act.ranks.size() || i < exp.tasks.size() || i < exp.ranks.size();
         ++i) {
      std::cerr << i << "\t";
      if (i < taskSizes.size()) {
        std::cerr << taskSizes[i] << "\t";
      } else {
        std::cerr << " "
                  << "\t";
      }

      if (i < exp.tasks.size()) {
        std::cerr << exp.tasks[i] << "\t";
      } else {
        std::cerr << " "
                  << "\t";
      }

      if (i < act.tasks.size()) {
        std::cerr << act.tasks[i] << "\t";
      } else {
        std::cerr << " "
                  << "\t";
      }

      if (i < exp.ranks.size()) {
        std::cerr << exp.ranks[i] << "\t";
      } else {
        std::cerr << " "
                  << "\t";
      }

      if (i < exp.ranks.size()) {
        std::cerr << exp.ranks[i] << "\t";
      } else {
        std::cerr << " "
                  << "\t";
      }

      std::cerr << "\n";
    }
  }
  return 0 == err;
}

template <typename TasksType, typename RanksType, typename TaskSizesView>
struct ThreadLoadBalanceIntoFunctor {
  ThreadLoadBalanceIntoFunctor(const TasksType &tasks, const RanksType &ranks, const TaskSizesView &taskSizes)
      : tasks_(tasks), ranks_(ranks), taskSizes_(taskSizes) {}

  // one thread does the whole merge
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t t) const {
    if (0 == t) {
      KokkosGraph::Impl::load_balance_into_thread(tasks_, ranks_, taskSizes_);
    }
  }

  TasksType tasks_;
  RanksType ranks_;
  TaskSizesView taskSizes_;
};

/*! \brief only works on the host
 */
template <typename TaskSize, typename Device>
auto run_load_balance_into_thread(const std::vector<TaskSize> &taskSizes) {
  using execution_space     = typename Device::execution_space;
  using Policy              = Kokkos::RangePolicy<execution_space>;
  using tasksize_view_type  = Kokkos::View<TaskSize *, Device>;
  using TaskIndex           = typename tasksize_view_type::size_type;
  using taskindex_view_type = Kokkos::View<TaskIndex *, Device>;
  using const_tasksize_uview_type =
      Kokkos::View<const TaskSize *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using tasksize_uview_type  = Kokkos::View<TaskSize *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using taskindex_uview_type = Kokkos::View<TaskIndex *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  TaskSize totalWork = 0;
  for (const auto &size : taskSizes) {
    totalWork += size;
  }

  // create device views of input data
  const_tasksize_uview_type uts(taskSizes.data(), taskSizes.size());
  tasksize_view_type ts("ts", uts.size());
  Kokkos::deep_copy(ts, uts);

  // device view of output data
  taskindex_view_type tasks("tasks", totalWork);
  tasksize_view_type ranks("ranks", totalWork);

  // run
  Kokkos::parallel_for(Policy(0, 1) /*one thread will run*/, ThreadLoadBalanceIntoFunctor(tasks, ranks, ts));

  // return output to host
  LBR<TaskSize, TaskIndex> lbr;
  lbr.tasks.resize(tasks.size());
  lbr.ranks.resize(ranks.size());
  taskindex_uview_type utasks(lbr.tasks.data(), lbr.tasks.size());
  tasksize_uview_type uranks(lbr.ranks.data(), lbr.ranks.size());
  Kokkos::deep_copy(utasks, tasks);
  Kokkos::deep_copy(uranks, ranks);
  Kokkos::fence();

  return lbr;
}

template <typename TasksType, typename RanksType, typename TaskSizesView>
struct TeamLoadBalanceIntoFunctor {
  TeamLoadBalanceIntoFunctor(const TasksType &tasks, const RanksType &ranks, const TaskSizesView &scratch,
                             const TaskSizesView &taskSizes)
      : tasks_(tasks), ranks_(ranks), scratch_(scratch), taskSizes_(taskSizes) {}

  // one team does the whole merge
  template <typename TeamMember>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember &handle) const {
    inclusive_prefix_sum_team(handle, scratch_, taskSizes_);
    KokkosGraph::load_balance_into_team(handle, tasks_, ranks_, scratch_);
  }

  TasksType tasks_;
  RanksType ranks_;
  TaskSizesView scratch_;
  TaskSizesView taskSizes_;
};

/*! \brief only works on the host
 */
template <typename TaskSize, typename Device>
auto run_load_balance_into_team(const std::vector<TaskSize> &taskSizes) {
  using execution_space     = typename Device::execution_space;
  using Policy              = Kokkos::TeamPolicy<execution_space>;
  using tasksize_view_type  = Kokkos::View<TaskSize *, Device>;
  using TaskIndex           = typename tasksize_view_type::size_type;
  using taskindex_view_type = Kokkos::View<TaskIndex *, Device>;
  using const_tasksize_uview_type =
      Kokkos::View<const TaskSize *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using tasksize_uview_type  = Kokkos::View<TaskSize *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using taskindex_uview_type = Kokkos::View<TaskIndex *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  TaskSize totalWork = 0;
  for (const auto &size : taskSizes) {
    totalWork += size;
  }

  // create device views of input data
  const_tasksize_uview_type uts(taskSizes.data(), taskSizes.size());
  tasksize_view_type ts("ts", uts.size());
  Kokkos::deep_copy(ts, uts);

  // device view of scratch data
  tasksize_view_type scratch("scratch", ts.size());

  // device view of output data
  taskindex_view_type tasks("tasks", totalWork);
  tasksize_view_type ranks("ranks", totalWork);

  // one team, small for non-GPU spaces
  const int leagueSize = 1;
  int teamSize;
  if constexpr (KokkosKernels::Impl::is_gpu_exec_space_v<execution_space>) {
    teamSize = 64;
  } else {
    teamSize = 1;
  }

  // run
  Kokkos::parallel_for(Policy(leagueSize, teamSize), TeamLoadBalanceIntoFunctor(tasks, ranks, scratch, ts));

  // return output to host
  LBR<TaskSize, TaskIndex> lbr;
  lbr.tasks.resize(tasks.size());
  lbr.ranks.resize(ranks.size());
  taskindex_uview_type utasks(lbr.tasks.data(), lbr.tasks.size());
  tasksize_uview_type uranks(lbr.ranks.data(), lbr.ranks.size());
  Kokkos::deep_copy(utasks, tasks);
  Kokkos::deep_copy(uranks, ranks);
  Kokkos::fence();

  return lbr;
}

/*! \brief only works on the host
 */
template <typename TaskSize, typename Device>
auto run_global_load_balance(const std::vector<TaskSize> &taskSizes) {
  using tasksize_view_type = Kokkos::View<TaskSize *, Device>;
  using TaskIndex          = typename tasksize_view_type::size_type;
  using const_tasksize_uview_type =
      Kokkos::View<const TaskSize *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using tasksize_uview_type  = Kokkos::View<TaskSize *, Device, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using taskindex_uview_type = Kokkos::View<TaskIndex *, Device, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // create device views of input data
  const_tasksize_uview_type uts(taskSizes.data(), taskSizes.size());
  tasksize_view_type ts("ts", uts.size());
  Kokkos::deep_copy(ts, uts);

  // run
  KokkosGraph::LoadBalanceResult<tasksize_view_type> lbr = KokkosGraph::load_balance(ts);
  Kokkos::fence();

  // return output to host
  LBR<TaskSize, TaskIndex> _lbr;
  _lbr.tasks.resize(lbr.tasks.size());
  _lbr.ranks.resize(lbr.ranks.size());
  taskindex_uview_type utasks(_lbr.tasks.data(), _lbr.tasks.size());
  tasksize_uview_type uranks(_lbr.ranks.data(), _lbr.ranks.size());
  Kokkos::deep_copy(utasks, lbr.tasks);
  Kokkos::deep_copy(uranks, lbr.ranks);
  Kokkos::fence();

  return _lbr;
}

template <typename TaskSize, typename Device>
void test_load_balance_into_thread(const std::vector<TaskSize> &taskSizes) {
  using TaskIndex = typename Kokkos::View<TaskSize *, Device>::size_type;

  // expected values
  LBR<TaskSize, TaskIndex> exp = vector_load_balance<TaskSize, TaskIndex>(taskSizes);

  // result under test
  LBR<TaskSize, TaskIndex> act = run_load_balance_into_thread<TaskSize, Device>(taskSizes);

  EXPECT_TRUE(vector_compare(taskSizes, exp, act));
}

template <typename TaskSize, typename Device>
void test_load_balance_into_team(const std::vector<TaskSize> &taskSizes) {
  using TaskIndex = typename Kokkos::View<TaskSize *, Device>::size_type;

  // expected values
  LBR<TaskSize, TaskIndex> exp = vector_load_balance<TaskSize, TaskIndex>(taskSizes);

  // result under test
  LBR<TaskSize, TaskIndex> act = run_load_balance_into_team<TaskSize, Device>(taskSizes);

  EXPECT_TRUE(vector_compare(taskSizes, exp, act));
}

template <typename TaskSize, typename Device>
void test_global_load_balance(const std::vector<TaskSize> &taskSizes) {
  using TaskIndex = typename Kokkos::View<TaskSize *, Device>::size_type;

  // expected values
  LBR<TaskSize, TaskIndex> exp = vector_load_balance<TaskSize, TaskIndex>(taskSizes);

  // result under test
  LBR<TaskSize, TaskIndex> act = run_global_load_balance<TaskSize, Device>(taskSizes);

  EXPECT_TRUE(vector_compare(taskSizes, exp, act));
}

/* test merge on the device, and compare to a simple host implementation
 */
template <typename TaskSize, typename Device>
void test_load_balance(const std::vector<TaskSize> &taskSizes) {
  test_load_balance_into_thread<TaskSize, Device>(taskSizes);
  test_load_balance_into_team<TaskSize, Device>(taskSizes);
  test_global_load_balance<TaskSize, Device>(taskSizes);
}

/* define specific and random merge test cases
 */
template <typename TaskSize, typename Device>
void test_load_balance() {
  test_load_balance<TaskSize, Device>({});
  test_load_balance<TaskSize, Device>({0});
  test_load_balance<TaskSize, Device>({1});
  test_load_balance<TaskSize, Device>({2});
  test_load_balance<TaskSize, Device>({0, 0});
  test_load_balance<TaskSize, Device>({0, 1});
  test_load_balance<TaskSize, Device>({1, 0});
  test_load_balance<TaskSize, Device>({0, 2});
  test_load_balance<TaskSize, Device>({2, 0});
  test_load_balance<TaskSize, Device>({3, 0, 2});
  test_load_balance<TaskSize, Device>({2, 3, 0});
  test_load_balance<TaskSize, Device>({2, 0, 2, 4});
  test_load_balance<TaskSize, Device>({7, 0, 13, 1, 11, 0});
  test_load_balance<TaskSize, Device>({0, 0, 0, 0, 1});

  for (size_t sz : {10, 100, 1000}) {
    std::vector<TaskSize> wi;
    for (size_t i = 0; i < sz; ++i) {
      wi.push_back(rand() % 20);
    }
    test_load_balance<TaskSize, Device>(wi);
  }
}

#define EXECUTE_TEST(TASK_SIZE, TASK_INDEX, DEVICE) \
  TEST_F(TestCategory, graph##_##load_balance##_##TASK_SIZE##_##DEVICE) { test_load_balance<TASK_SIZE, DEVICE>(); }

// TODO tests with instantiated offset types

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(int, int, TestDevice)
#endif

#if 0
#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(int64_t, int64_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(size_t, size_t, TestDevice)
#endif
#endif

#undef EXECUTE_TEST
