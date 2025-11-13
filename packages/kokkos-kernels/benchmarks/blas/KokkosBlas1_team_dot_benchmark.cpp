// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>

#include <Kokkos_Core.hpp>
#include <KokkosBlas1_team_dot.hpp>
#include <Kokkos_Random.hpp>
#include "KokkosKernels_TestStringUtils.hpp"

#include <benchmark/benchmark.h>

// Functor to handle the case of a "without Cuda" build
template <class Vector, class ExecSpace>
struct teamDotFunctor {
  // Compile - time check to see if your data type is a Kokkos::View:
  static_assert(Kokkos::is_view<Vector>::value, "Vector is not a Kokkos::View.");

  using Scalar = typename Vector::non_const_value_type;
  // Vector is templated on memory space
  using execution_space = ExecSpace;  // Kokkos Execution Space
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;

  // Declare Kokkos::View Vectors, x and y
  Vector x;
  Vector y;

  // Functor instead of KOKKOS_LAMBDA expression

  KOKKOS_INLINE_FUNCTION void operator()(const team_member& team) const { KokkosBlas::Experimental::dot(team, x, y); }
  // Constructor
  teamDotFunctor(Vector X_, Vector Y_) {
    x = X_;
    y = Y_;
  }
};

template <class ExecSpace>
static void run(benchmark::State& state) {
  const auto m = state.range(0);
  // Declare type aliases
  using Scalar   = double;
  using MemSpace = typename ExecSpace::memory_space;

  // For the Team implementation of dot; ExecSpace is implicit;
  using policy = Kokkos::TeamPolicy<ExecSpace>;

  // Create 1D view w/ Device as the ExecSpace; this is an input vector
  Kokkos::View<Scalar*, MemSpace> x("X", m);
  // Create 1D view w/ Device as the ExecSpace; this is the output vector
  Kokkos::View<Scalar*, MemSpace> y("Y", m);

  // Here, deep_copy is filling / copying values into Host memory from Views X
  // and Y
  Kokkos::deep_copy(x, 3.0);
  Kokkos::deep_copy(y, 2.0);

  Kokkos::fence();
  for (auto _ : state) {
    teamDotFunctor<Kokkos::View<Scalar*, MemSpace>, ExecSpace> teamDotFunctorLiveTestInstance(x, y);
    Kokkos::parallel_for("TeamDotUsage -- Live Test", policy(1, Kokkos::AUTO), teamDotFunctorLiveTestInstance);
    Kokkos::fence();
  }
  const size_t iterFlop   = (size_t)2 * m;
  const size_t totalFlop  = iterFlop * state.iterations();
  state.counters["FLOP"]  = benchmark::Counter(iterFlop);
  state.counters["FLOPS"] = benchmark::Counter(totalFlop, benchmark::Counter::kIsRate);
}

BENCHMARK(run<Kokkos::DefaultExecutionSpace>)
    ->Name("KokkosBlas_team_dot/run<Kokkos::DefaultExecutionSpace>")
    ->Unit(benchmark::kMicrosecond)
    ->UseRealTime()
    ->ArgName("m")
    ->RangeMultiplier(10)
    ->Range(100000, 100000000);
