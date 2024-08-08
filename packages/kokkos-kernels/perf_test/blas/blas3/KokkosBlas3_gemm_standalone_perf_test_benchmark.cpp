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

#include "KokkosBlas3_gemm.hpp"
#include <Kokkos_Random.hpp>
#include "KokkosKernels_TestUtils.hpp"
#include "KokkosKernels_perf_test_utilities.hpp"
#include "Benchmark_Context.hpp"
#include <benchmark/benchmark.h>

struct blas3_gemm_params : public perf_test::CommonInputParams {
  int m = 1000;
  int n = 1000;
  int k = 1000;

  static blas3_gemm_params get_params(int& argc, char** argv) {
    blas3_gemm_params params;
    perf_test::parse_common_options(argc, argv, params);

    for (int i = 1; i < argc; ++i) {
      if (perf_test::check_arg_int(i, argc, argv, "--m", params.m)) {
        ++i;
      } else if (perf_test::check_arg_int(i, argc, argv, "--n", params.n)) {
        ++i;
      } else if (perf_test::check_arg_int(i, argc, argv, "--k", params.k)) {
        ++i;
      } else if (std::string(argv[i]).find("--benchmark") == 0) {
        continue;  // ignore benchmark arguments
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

    std::cerr << "\t[Optional] --m      :: Rows in A (default 1000)" << std::endl;
    std::cerr << "\t[Optional] --n      :: Columns in A / Rows in B (default 1000)" << std::endl;
    std::cerr << "\t[Optional] --k      :: Columns in B (default 1000)" << std::endl;
  }
};

template <typename ExecSpace, typename Scalar, typename ALayout, typename BLayout>
static void KokkosBlas3_GEMM(benchmark::State& state) {
  const auto m = state.range(0);
  const auto n = state.range(1);
  const auto k = state.range(2);

  using MemSpace = typename ExecSpace::memory_space;
  using Device   = Kokkos::Device<ExecSpace, MemSpace>;
  Kokkos::View<Scalar**, ALayout, Device> A(Kokkos::view_alloc(Kokkos::WithoutInitializing, "A"), m, n);
  Kokkos::View<Scalar**, BLayout, Device> B(Kokkos::view_alloc(Kokkos::WithoutInitializing, "B"), n, k);
  Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device> C(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C"), m, k);
  Kokkos::Random_XorShift64_Pool<ExecSpace> pool(123);
  Kokkos::fill_random(A, pool, 10.0);
  Kokkos::fill_random(B, pool, 10.0);

  // Do a warm-up run
  KokkosBlas::gemm("N", "N", 1.0, A, B, 0.0, C);
  Kokkos::fence();
  double total_time = 0.0;

  for (auto _ : state) {
    Kokkos::Timer timer;
    KokkosBlas::gemm("N", "N", 1.0, A, B, 0.0, C);
    ExecSpace().fence();

    double time = timer.seconds();
    total_time += time;
    state.SetIterationTime(time);
  }

  state.counters[ExecSpace::name()]    = 1;
  state.counters["Avg GEMM time (s):"] = benchmark::Counter(total_time, benchmark::Counter::kAvgIterations);
  size_t flopsPerRun                   = (size_t)2 * m * n * k;
  state.counters["Avg GEMM (FLOP/s):"] = benchmark::Counter(flopsPerRun, benchmark::Counter::kIsIterationInvariantRate);
  if constexpr (std::is_same_v<ALayout, Kokkos::LayoutLeft>) {
    state.counters["Memory Layout in A: LayoutLeft"] = 1;
  } else {
    state.counters["Memory Layout in A: LayoutRight"] = 1;
  }
  if constexpr (std::is_same_v<BLayout, Kokkos::LayoutLeft>) {
    state.counters["Memory Layout in B: LayoutLeft"] = 1;
  } else {
    state.counters["Memory Layout in B: LayoutRight"] = 1;
  }
}

template <typename ExecSpace>
void run(const blas3_gemm_params& params) {
  using LL     = Kokkos::LayoutLeft;
  using LR     = Kokkos::LayoutRight;
  using Scalar = double;

  const auto name      = "KokkosBlas3_GEMM";
  const auto arg_names = std::vector<std::string>{"m", "n", "k"};
  const auto args      = std::vector<int64_t>{params.m, params.n, params.k};

  KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas3_GEMM<ExecSpace, Scalar, LL, LL>, arg_names, args,
                                             params.repeat);
  KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas3_GEMM<ExecSpace, Scalar, LL, LR>, arg_names, args,
                                             params.repeat);
  KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas3_GEMM<ExecSpace, Scalar, LR, LL>, arg_names, args,
                                             params.repeat);
  KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas3_GEMM<ExecSpace, Scalar, LR, LR>, arg_names, args,
                                             params.repeat);
}

int main(int argc, char** argv) {
  const auto params     = blas3_gemm_params::get_params(argc, argv);
  const int num_threads = params.use_openmp;

  // the common parameter parser takes the requested device ID and
  // adds 1 to it (e.g. --cuda 0 -> params.use_cuda = 1)
  // this is presumably so that 0 can be a sentinel value,
  // even though device ID 0 is valid
  // here, we use CUDA, SYCL, or HIP, whichever is set first, to
  // choose which device Kokkos should initialize on
  // or -1, for no such selection
  const int device_id = params.use_cuda
                            ? params.use_cuda - 1
                            : (params.use_sycl ? params.use_sycl - 1 : (params.use_hip ? params.use_hip - 1 : -1));

  Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(num_threads).set_device_id(device_id));
  benchmark::Initialize(&argc, argv);
  benchmark::SetDefaultTimeUnit(benchmark::kSecond);
  KokkosKernelsBenchmark::add_benchmark_context(true);

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
