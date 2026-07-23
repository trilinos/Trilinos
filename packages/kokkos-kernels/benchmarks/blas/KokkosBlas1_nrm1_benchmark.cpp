// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include "KokkosBlas1_nrm1.hpp"
#include <benchmark/benchmark.h>

template <class ExecSpace>
static void run(benchmark::State& state) {
  const auto m = state.range(0);
  // Declare type aliases
  using Scalar   = double;
  using MemSpace = typename ExecSpace::memory_space;
  using Device   = Kokkos::Device<ExecSpace, MemSpace>;

  // Create 1D view w/ Device as the ExecSpace; this is an input vector
  // A(view_alloc(WithoutInitializing, "label"), m, n);
  Kokkos::View<Scalar*, Device> x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "x"), m);

  // Create 1D view w/ Device as the ExecSpace; this is the output vector
  Kokkos::View<Scalar, Device> r("result");

  // Declaring variable pool w/ a seeded random number;
  // a parallel random number generator, so you
  // won't get the same number with a given seed each time
  Kokkos::Random_XorShift64_Pool<ExecSpace> pool(123);

  Kokkos::fill_random(x, pool, 10.0);
  ExecSpace space{};

  Kokkos::fence();
  for (auto _ : state) {
    KokkosBlas::nrm1(space, r, x);
    space.fence();
  }
  const size_t iterFlop   = (size_t)2 * m;
  const size_t totalFlop  = iterFlop * state.iterations();
  state.counters["FLOP"]  = benchmark::Counter(iterFlop);
  state.counters["FLOPS"] = benchmark::Counter(totalFlop, benchmark::Counter::kIsRate);
}

BENCHMARK(run<Kokkos::DefaultExecutionSpace>)
    ->Name("KokkosBlas_nrm1")
    ->Unit(benchmark::kMicrosecond)
    ->UseRealTime()
    ->ArgName("m")
    ->RangeMultiplier(10)
    ->Range(100000, 100000000);
