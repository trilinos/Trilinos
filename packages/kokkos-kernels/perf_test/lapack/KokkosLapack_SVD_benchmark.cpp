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

#include "KokkosLapack_svd.hpp"

#include "KokkosKernels_TestUtils.hpp"
#include "KokkosKernels_perf_test_utilities.hpp"

#include <benchmark/benchmark.h>
#include "Benchmark_Context.hpp"

struct svd_parameters {
  int numRows, numCols;
  bool verbose;

  svd_parameters(const int numRows_, const int numCols_, const bool verbose_)
      : numRows(numRows_), numCols(numCols_), verbose(verbose_){};
};

void print_options() {
  std::cerr << "Options\n" << std::endl;

  std::cerr << perf_test::list_common_options();

  std::cerr << "\t[Optional] --verbose     :: enable verbose output" << std::endl;
  std::cerr << "\t[Optional] --m           :: number of rows of A" << std::endl;
  std::cerr << "\t[Optional] --n           :: number of columns of A" << std::endl;
}  // print_options

int parse_inputs(svd_parameters& params, int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    if (perf_test::check_arg_int(i, argc, argv, "--m", params.numRows)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--n", params.numCols)) {
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--verbose", params.verbose)) {
    } else {
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      print_options();
      return 1;
    }
  }
  return 0;
}  // parse_inputs

template <class ExecutionSpace>
void run_svd_benchmark(benchmark::State& state, const svd_parameters& svd_params) {
  using mat_type = Kokkos::View<double**, Kokkos::LayoutLeft, ExecutionSpace>;
  using vec_type = Kokkos::View<double*, Kokkos::LayoutLeft, ExecutionSpace>;

  const int m = svd_params.numRows;
  const int n = svd_params.numCols;

  mat_type A("A", m, n), U("U", m, m), Vt("Vt", n, n);
  vec_type S("S", Kokkos::min(m, n));

  const uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  Kokkos::Random_XorShift64_Pool<ExecutionSpace> rand_pool(seed);

  // Initialize A with random numbers
  double randStart = 0, randEnd = 0;
  Test::getRandomBounds(10.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);

  for (auto _ : state) {
    (void)_;
    KokkosLapack::svd("A", "A", A, S, U, Vt);
    Kokkos::fence();
  }
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);

  benchmark::Initialize(&argc, argv);
  benchmark::SetDefaultTimeUnit(benchmark::kMillisecond);
  KokkosKernelsBenchmark::add_benchmark_context(true);

  perf_test::CommonInputParams common_params;
  perf_test::parse_common_options(argc, argv, common_params);
  svd_parameters svd_params(0, 0, false);
  parse_inputs(svd_params, argc, argv);

  std::string bench_name = "KokkosLapack_SVD";

  if (0 < common_params.repeat) {
    benchmark::RegisterBenchmark(bench_name.c_str(), run_svd_benchmark<Kokkos::DefaultExecutionSpace>, svd_params)
        ->UseRealTime()
        ->Iterations(common_params.repeat);
  } else {
    benchmark::RegisterBenchmark(bench_name.c_str(), run_svd_benchmark<Kokkos::DefaultExecutionSpace>, svd_params)
        ->UseRealTime();
  }

  benchmark::RunSpecifiedBenchmarks();

  benchmark::Shutdown();
  Kokkos::finalize();

  return 0;
}
