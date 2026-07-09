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

#include "KokkosKernels_ExecSpaceUtils.hpp"

#include "KokkosGraph_MergePath.hpp"

struct StepperParams {
  KokkosGraph::StepDirection dir;
  KokkosGraph::StepperContext ctx;
};

// lexical compare of StepperContext for std::set
KOKKOS_INLINE_FUNCTION
bool operator<(const KokkosGraph::StepperContext &a, const KokkosGraph::StepperContext &b) {
  if (a.ai < b.ai) {
    return true;
  } else if (a.ai > b.ai) {
    return false;
  } else {
    if (a.bi < b.bi) {
      return true;
    } else if (a.bi > b.bi) {
      return false;
    } else {
      return a.pi < b.pi;
    }
  }
}

// lexical compare of StepperParams for std::set
KOKKOS_INLINE_FUNCTION
bool operator<(const StepperParams &a, const StepperParams &b) {
  if (a.dir < b.dir) {
    return true;
  } else if (a.dir > b.dir) {
    return false;
  } else {
    return a.ctx < b.ctx;
  }
}

// store all parameters provided to the stepper
template <typename AView, typename BView, typename ParamView>
struct Stepper {
  KOKKOS_INLINE_FUNCTION
  Stepper(const AView &_a, const BView &_b, const ParamView &_params) : a(_a), b(_b), params(_params) {}

  /// \brief stepper for thread
  KOKKOS_INLINE_FUNCTION
  void operator()(const KokkosGraph::StepDirection &dir, const KokkosGraph::StepperContext &step) const {
    params(step.pi) = {dir, step};
  }

  /// \brief stepper for team
  KOKKOS_INLINE_FUNCTION
  void operator()(const KokkosGraph::StepDirection &dir, const KokkosGraph::StepperContext &step,
                  const KokkosGraph::StepperContext &thread) const {
    KokkosGraph::StepperContext total(step.ai + thread.ai, step.bi + thread.bi, step.pi + thread.pi);
    params(thread.pi + step.pi) = {dir, total};
  }

  AView a;
  BView b;
  ParamView params;
};

/*! \brief Calls the Stepper in a team context
 */
template <typename AView, typename BView, typename ParamView>
struct TeamCaller {
  TeamCaller(const AView &_a, const BView &_b, const ParamView &_params) : a(_a), b(_b), params(_params) {}

  template <typename Member>
  KOKKOS_INLINE_FUNCTION void operator()(const Member &handle) const {
    Stepper stepper(a, b, params);
    KokkosGraph::merge_path_team(handle, a, b, a.size() + b.size(), stepper);
  }

  AView a;
  BView b;
  ParamView params;
};

template <typename T, typename Device>
void test_merge_path_team(const std::vector<T> &_a, const std::vector<T> &_b,
                          const std::set<StepperParams> &expParams) {
  using execution_space = typename Device::execution_space;
  using Policy          = Kokkos::TeamPolicy<execution_space>;

  // copy input data to device
  Kokkos::View<const T *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> ua(_a.data(), _a.size());
  Kokkos::View<const T *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> ub(_b.data(), _b.size());
  Kokkos::View<T *, Device> a("", ua.size());
  Kokkos::View<T *, Device> b("", ub.size());
  Kokkos::deep_copy(a, ua);
  Kokkos::deep_copy(b, ub);

  // create actual params
  Kokkos::View<StepperParams *, Device> actParams("", a.size() + b.size());

  // one team, small for non-GPU spaces
  const int leagueSize = 1;
  int teamSize;
  Policy policy;
  if constexpr (KokkosKernels::Impl::is_gpu_exec_space_v<execution_space>) {
    teamSize = 64;
  } else {
    teamSize = 1;
  }

  // run merge
  Kokkos::parallel_for(Policy(leagueSize, teamSize), TeamCaller(a, b, actParams));
  Kokkos::fence();

  // copy results back to host
  Kokkos::View<StepperParams *, Kokkos::HostSpace> hostParams("", actParams.size());
  Kokkos::deep_copy(hostParams, actParams);

  // every parameter set visited once
  auto remParams = expParams;
  for (size_t i = 0; i < hostParams.size(); ++i) {
    StepperParams needle = hostParams(i);
    EXPECT_TRUE(remParams.count(needle) > 0);
    remParams.erase(needle);
  }
  // nothing unvisited
  EXPECT_TRUE(remParams.empty());
}

/*! \brief Calls the Stepper in a thread context
 */
template <typename AView, typename BView, typename ParamView>
struct ThreadCaller {
  ThreadCaller(const AView &_a, const BView &_b, const ParamView &_params) : a(_a), b(_b), params(_params) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (0 == i) {
      Stepper stepper(a, b, params);
      KokkosGraph::merge_path_thread(a, b, a.size() + b.size(), stepper);
    }
  }

  AView a;
  BView b;
  ParamView params;
};

template <typename T, typename Device>
void test_merge_path_thread(const std::vector<T> &_a, const std::vector<T> &_b, std::set<StepperParams> expParams) {
  using execution_space = typename Device::execution_space;
  using Policy          = Kokkos::RangePolicy<execution_space>;

  // copy input data to device
  Kokkos::View<const T *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> ua(_a.data(), _a.size());
  Kokkos::View<const T *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> ub(_b.data(), _b.size());
  Kokkos::View<T *, Device> a("", ua.size());
  Kokkos::View<T *, Device> b("", ub.size());
  Kokkos::deep_copy(a, ua);
  Kokkos::deep_copy(b, ub);

  // place to record params passed to the stepper function
  Kokkos::View<StepperParams *, Device> actParams("", a.size() + b.size());

  // run
  Kokkos::parallel_for(Policy(0, 1), ThreadCaller(a, b, actParams));
  Kokkos::fence();

  // copy results back to host
  Kokkos::View<StepperParams *, Kokkos::HostSpace> hostParams("", actParams.size());
  Kokkos::deep_copy(hostParams, actParams);

  // every parameter set visited once
  for (size_t i = 0; i < hostParams.size(); ++i) {
    StepperParams needle = hostParams(i);
    EXPECT_TRUE(expParams.count(needle) > 0);
    expParams.erase(needle);
  }
  // nothing unvisited
  EXPECT_TRUE(expParams.empty());
}

template <typename T, typename Device>
void test_merge_path(const std::vector<T> &a, const std::vector<T> &b, const std::set<StepperParams> &expParams) {
  test_merge_path_thread<T, Device>(a, b, expParams);
  test_merge_path_team<T, Device>(a, b, expParams);
}

template <typename T, typename Device>
void test_merge_path() {
  constexpr KokkosGraph::StepDirection B = KokkosGraph::StepDirection::b;
  constexpr KokkosGraph::StepDirection A = KokkosGraph::StepDirection::a;
  (void)A;  // warning: variable "A" was unused, g++ 10.1
  (void)B;  // warning: variable "B" was unused, g++ 10.1

  // empty lists
  test_merge_path<T, Device>({},  // a
                             {},  // b
                             {}   // expected parameters to stepper function
  );

  // B
  test_merge_path<T, Device>({}, {0}, {{B, {0, 0, 0}}});

  // A
  test_merge_path<T, Device>({0}, {}, {{A, {0, 0, 0}}});

  // A, B
  test_merge_path<T, Device>({0}, {0}, {{A, {0, 0, 0}}, {B, {1, 0, 1}}});

  // A, B, A
  test_merge_path<T, Device>({0, 1}, {0}, {{A, {0, 0, 0}}, {B, {1, 0, 1}}, {A, {1, 1, 2}}});

  // A, A, B, B
  test_merge_path<T, Device>({0, 0}, {0, 0}, {{A, {0, 0, 0}}, {A, {1, 0, 1}}, {B, {2, 0, 2}}, {B, {2, 1, 3}}});
}

#define EXECUTE_TEST(ORDINAL, DEVICE) \
  TEST_F(TestCategory, graph##_##merge_path##_##ORDINAL##_##DEVICE) { test_merge_path<ORDINAL, DEVICE>(); }

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(int64_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(size_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_SCALAR_FLOAT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_SCALAR_DOUBLE)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, TestDevice)
#endif

#undef EXECUTE_TEST
