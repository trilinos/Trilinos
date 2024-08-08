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

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include "KokkosBlas2_gemv.hpp"

#include "KokkosKernels_TestUtils.hpp"
#include "KokkosKernels_perf_test_utilities.hpp"

#include <Benchmark_Context.hpp>
#include <benchmark/benchmark.h>

struct blas2_gemv_params : public perf_test::CommonInputParams {
  int m           = 5000;
  int n           = 5000;
  bool layoutLeft = true;

  static blas2_gemv_params get_params(int& argc, char** argv) {
    blas2_gemv_params params;
    perf_test::parse_common_options(argc, argv, params);

    for (int i = 1; i < argc; ++i) {
      if (perf_test::check_arg_int(i, argc, argv, "--m", params.m)) {
        ++i;
      } else if (perf_test::check_arg_int(i, argc, argv, "--n", params.n)) {
        ++i;
      } else if (std::string layout; perf_test::check_arg_str(i, argc, argv, "--layout", layout)) {
        if (0 == Test::string_compare_no_case(layout, "left"))
          params.layoutLeft = true;
        else if (0 == Test::string_compare_no_case(layout, "right"))
          params.layoutLeft = false;
        else {
          std::cerr << "Invalid layout: must be 'left' or 'right'.\n";
          exit(1);
        }
        ++i;
      } else {
        std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
        print_options();
        exit(1);
      }
    }
    return params;
  }

  static void print_options() {
    std::cerr << "Options\n" << std::endl;
    std::cerr << perf_test::list_common_options();

    std::cerr << "\t[Optional] --m      :: number of rows to generate (default 5000)" << std::endl;
    std::cerr << "\t[Optional] --n      :: number of cols to generate (default 5000)" << std::endl;
    std::cerr << "\t[Optional] --layout :: matrix layout ('left' or 'right', "
                 "default 'left')"
              << std::endl;
  }
};

template <typename Scalar, typename Layout, typename ExecSpace>
static void KokkosBlas2_GEMV(benchmark::State& state) {
  const auto m = state.range(0);
  const auto n = state.range(1);

  // Declare type aliases
  using MemSpace = typename ExecSpace::memory_space;
  using Device   = Kokkos::Device<ExecSpace, MemSpace>;

  // Create a View containing a 2D matrix; allocate KokkosView with template
  // args of Scalar**, a layout, and
  Kokkos::View<Scalar**, Layout, Device> A(Kokkos::view_alloc(Kokkos::WithoutInitializing, "A"), m, n);
  // Create Views containing 1D matrix; allocate (without) matrix "x" of size n
  Kokkos::View<Scalar*, Device> x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "x"), n);
  // Create Views containing 1D matrix; allocate (without) matrix "y" of size m
  Kokkos::View<Scalar*, Device> y(Kokkos::view_alloc(Kokkos::WithoutInitializing, "y"), m);

  // Declaring variable pool w/ a number seed;
  // a parallel random number generator, so you
  // won't get the same number with a given seed each time
  Kokkos::Random_XorShift64_Pool<ExecSpace> pool(123);

  // Fill 2D Matrix "A" and 1D matrix (i.e., a vector) "x" with random values;
  // Here, 10 is the max value of the random generator between 1 and 10
  // (uniform )
  Kokkos::fill_random(A, pool, 10.0);
  Kokkos::fill_random(x, pool, 10.0);

  // Do a warm-up run
  KokkosBlas::gemv("N", 1.0, A, x, 0.0, y);
  Kokkos::fence();
  double total_time = 0.0;

  for (auto _ : state) {
    // Start timing
    Kokkos::Timer timer;
    KokkosBlas::gemv("N", 1.0, A, x, 0.0, y);
    ExecSpace().fence();

    double time = timer.seconds();
    total_time += time;
    state.SetIterationTime(time);
  }

  state.counters[ExecSpace::name()]    = 1;
  state.counters["Avg GEMV time (s):"] = benchmark::Counter(total_time, benchmark::Counter::kAvgIterations);
  size_t flopsPerRun                   = (size_t)2 * m * n;
  state.counters["Avg GEMV FLOP/s:"]   = benchmark::Counter(flopsPerRun, benchmark::Counter::kIsIterationInvariantRate);
}

template <typename ExecSpace>
void run(const blas2_gemv_params& params) {
  using Scalar = double;

  const auto name      = "KokkosBlas2_GEMV";
  const auto arg_names = std::vector<std::string>{"m", "n", params.layoutLeft ? "LayoutLeft" : "LayoutRight"};
  const auto args      = std::vector<int64_t>{params.m, params.n, 1};

  if (params.layoutLeft) {
    KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas2_GEMV<Scalar, Kokkos::LayoutLeft, ExecSpace>, arg_names,
                                               args, params.repeat);
  } else {
    KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas2_GEMV<Scalar, Kokkos::LayoutRight, ExecSpace>,
                                               arg_names, args, params.repeat);
  }
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  benchmark::Initialize(&argc, argv);
  benchmark::SetDefaultTimeUnit(benchmark::kSecond);
  KokkosKernelsBenchmark::add_benchmark_context(true);

  const auto params = blas2_gemv_params::get_params(argc, argv);

  if (params.use_threads) {
#if defined(KOKKOS_ENABLE_THREADS)
    run<Kokkos::Threads>(params);
#else
    std::cout << "ERROR:  PThreads requested, but not available.\n";
    return 1;
#endif
  }

  if (params.use_openmp) {
#if defined(KOKKOS_ENABLE_OPENMP)
    run<Kokkos::OpenMP>(params);
#else
    std::cout << "ERROR: OpenMP requested, but not available.\n";
    return 1;
#endif
  }

  if (params.use_cuda) {
#if defined(KOKKOS_ENABLE_CUDA)
    run<Kokkos::Cuda>(params);
#else
    std::cout << "ERROR: CUDA requested, but not available.\n";
    return 1;
#endif
  }

  if (params.use_hip) {
#if defined(KOKKOS_ENABLE_HIP)
    run<Kokkos::HIP>(params);
#else
    std::cout << "ERROR: HIP requested, but not available.\n";
    return 1;
#endif
  }

  if (params.use_sycl) {
#if defined(KOKKOS_ENABLE_SYCL)
    run<Kokkos::Experimental::SYCL>(params);
#else
    std::cout << "ERROR: SYCL requested, but not available.\n";
    return 1;
#endif
  }

  // use serial if no backend is specified
  if (!params.use_cuda and !params.use_hip and !params.use_openmp and !params.use_sycl and !params.use_threads) {
#if defined(KOKKOS_ENABLE_SERIAL)
    run<Kokkos::Serial>(params);
#else
    std::cout << "ERROR: Serial device requested, but not available.\n";
    return 1;
#endif
  }

  benchmark::RunSpecifiedBenchmarks();

  benchmark::Shutdown();
  Kokkos::finalize();
  return 0;
}
