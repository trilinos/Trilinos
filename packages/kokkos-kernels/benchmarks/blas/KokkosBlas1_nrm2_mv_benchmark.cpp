// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <KokkosBlas1_nrm2.hpp>
#include <benchmark/benchmark.h>

template <class ExecSpace>
static void run(benchmark::State& state) {
  const auto m = state.range(0);
  const auto n = state.range(1);
  // Declare type aliases
  using Scalar   = double;
  using MemSpace = typename ExecSpace::memory_space;
  using Device   = Kokkos::Device<ExecSpace, MemSpace>;

  Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device> x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "x"), m, n);

  Kokkos::View<Scalar*, Device> norms("norms", n);

  // Declaring variable pool w/ a seeded random number;
  // a parallel random number generator, so you
  // won't get the same number with a given seed each time
  Kokkos::Random_XorShift64_Pool<ExecSpace> pool(123);

  Kokkos::fill_random(x, pool, 10.0);
  ExecSpace space;

  Kokkos::fence();
  for (auto _ : state) {
    KokkosBlas::nrm2(space, norms, x);
    space.fence();
  }

  const size_t iterFlop   = (size_t)2 * m * n;
  const size_t totalFlop  = iterFlop * state.iterations();
  state.counters["FLOP"]  = benchmark::Counter(iterFlop);
  state.counters["FLOPS"] = benchmark::Counter(totalFlop, benchmark::Counter::kIsRate);
}

BENCHMARK(run<Kokkos::DefaultExecutionSpace>)
    ->Name("KokkosBlas_nrm2_mv")
    ->Unit(benchmark::kMicrosecond)
    ->UseRealTime()
    ->ArgNames({"m", "n"})
    ->ArgsProduct({benchmark::CreateRange(100000, 100000000, 10), benchmark::CreateRange(5, 5, 1)});
