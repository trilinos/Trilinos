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

#include "KokkosKernels_helpers.hpp"
#include "KokkosBlas2_ger.hpp"
#include <typeinfo>

#include "KokkosKernels_TestUtils.hpp"
#include "KokkosKernels_perf_test_utilities.hpp"

#include <Benchmark_Context.hpp>
#include <benchmark/benchmark.h>

struct blas2_ger_params : public perf_test::CommonInputParams {
  int verbosity          = 0;
  int m                  = 5000;
  int n                  = 5000;
  bool layoutLeft        = true;
  std::string scalarType = "double";
  std::string yMode      = "transpose";

  static blas2_ger_params get_params(int& argc, char** argv) {
    blas2_ger_params params;
    perf_test::parse_common_options(argc, argv, params);

    for (int i = 1; i < argc; ++i) {
      if (perf_test::check_arg_int(i, argc, argv, "--verbosity", params.verbosity)) {
        ++i;
      } else if (perf_test::check_arg_int(i, argc, argv, "--m", params.m)) {
        ++i;
      } else if (perf_test::check_arg_int(i, argc, argv, "--n", params.n)) {
        ++i;
      } else if (std::string layout; perf_test::check_arg_str(i, argc, argv, "--layout", layout)) {
        if (0 == Test::string_compare_no_case(layout, "left"))
          params.layoutLeft = true;
        else if (0 == Test::string_compare_no_case(layout, "right"))
          params.layoutLeft = false;
        else {
          std::cerr << "Invalid '--layout': must be 'left' or 'right'.\n";
          exit(1);
        }
        ++i;
      } else if (std::string scalarType; perf_test::check_arg_str(i, argc, argv, "--scalarType", scalarType)) {
        if ((0 == Test::string_compare_no_case(scalarType, "int32")) ||
            (0 == Test::string_compare_no_case(scalarType, "int64")) ||
            (0 == Test::string_compare_no_case(scalarType, "float")) ||
            (0 == Test::string_compare_no_case(scalarType, "double")) ||
            (0 == Test::string_compare_no_case(scalarType, "complex_float")) ||
            (0 == Test::string_compare_no_case(scalarType, "complex_double"))) {
          params.scalarType = scalarType;
        } else {
          std::cerr << "Invalid '--scalarType': must be 'int32' or 'int64' or "
                       "'float' or 'double' or 'complex_float' or "
                       "'complex_double'.\n";
          exit(1);
        }
        ++i;
      } else if (std::string yMode; perf_test::check_arg_str(i, argc, argv, "--yMode", yMode)) {
        if ((0 == Test::string_compare_no_case(yMode, "transpose")) ||
            (0 == Test::string_compare_no_case(yMode, "Hermitian"))) {
          params.yMode = yMode;
        } else {
          std::cerr << "Invalid '--yMode': must be 'transpose' or 'Hermitian'.\n";
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

    std::cerr << "\t[Optional] --m :: number of rows to generate (default 5000)" << std::endl;
    std::cerr << "\t[Optional] --n :: number of cols to generate (default 5000)" << std::endl;
    std::cerr << "\t[Optional] --layout :: matrix layout ('left' or 'right', "
                 "default 'left')"
              << std::endl;
    std::cerr << "\t[Optional] --scalarType :: scalar type ('int32' or 'int64'"
                 " or 'float' or 'double' or 'complex_float' or 'complex_double'"
                 ", default 'double')"
              << std::endl;
    std::cerr << "\t[Optional] --yMode :: y mode ('transpose' or 'Hermitian'"
                 ", default 'transpose')"
              << std::endl;
  }
};

template <typename tScalar, typename tLayout, typename tExecSpace>
static void KokkosBlas2_GER(benchmark::State& state) {
  const auto verbosity    = state.range(0);
  const auto m            = state.range(1);
  const auto n            = state.range(2);
  const auto yIsTranspose = state.range(3);
  tScalar a(Kokkos::ArithTraits<tScalar>::zero());

  if (verbosity > 0) {
    std::cout << "Entering KokkosBlas2_GER()"
              << ": m = " << m << ", n = " << n << ", yIsTranspose = " << yIsTranspose
              << ", tScalar = " << Kokkos::ArithTraits<tScalar>::name() << ", tLayout = " << typeid(tLayout).name()
              << std::endl;
  }

  using MemSpace = typename tExecSpace::memory_space;
  using Device   = Kokkos::Device<tExecSpace, MemSpace>;

  Kokkos::View<tScalar**, tLayout, Device> A(Kokkos::view_alloc(Kokkos::WithoutInitializing, "A"), m, n);
  Kokkos::View<tScalar*, Device> x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "x"), m);
  Kokkos::View<tScalar*, Device> y(Kokkos::view_alloc(Kokkos::WithoutInitializing, "y"), n);

  Kokkos::Random_XorShift64_Pool<tExecSpace> pool(123);

  char yMode('t');
  if (!yIsTranspose) yMode = 'H';

  tScalar rangeValue(Kokkos::ArithTraits<tScalar>::zero());
  if constexpr (Kokkos::ArithTraits<tScalar>::isOrdinal) {
    rangeValue = 10;
    a          = 3;
  } else if constexpr (Kokkos::ArithTraits<tScalar>::is_complex) {
    rangeValue.real() = 10.;
    rangeValue.imag() = 10.;
    a                 = tScalar(2.5, 3.6);
  } else {
    rangeValue = 10.;
    a          = 2.5;
  }

  // Fill 'A', 'x', and 'y' with samples from an uniform random variable with
  // range [1,rangeValue)
  Kokkos::fill_random(A, pool, rangeValue);
  Kokkos::fill_random(x, pool, rangeValue);
  Kokkos::fill_random(y, pool, rangeValue);

  if (verbosity > 0) {
    std::cout << "In KokkosBlas2_GER()"
              << ": yMode = " << yMode << ", a = " << a << ", rangeValue = " << rangeValue << std::endl;
  }

  // Do a warm-up run
  KokkosBlas::ger(&yMode, a, x, y, A);
  Kokkos::fence();
  double total_time = 0.0;

  for (auto _ : state) {
    // Start timing
    Kokkos::Timer timer;
    KokkosBlas::ger(&yMode, a, x, y, A);
    tExecSpace().fence();

    double time = timer.seconds();
    total_time += time;
    state.SetIterationTime(time);
  }

  state.counters[tExecSpace::name()]  = 1;
  state.counters["Avg GER time (s):"] = benchmark::Counter(total_time, benchmark::Counter::kAvgIterations);
  size_t flopsPerRun                  = (size_t)3 * m * n;
  state.counters["Avg GER FLOP/s:"]   = benchmark::Counter(flopsPerRun, benchmark::Counter::kIsIterationInvariantRate);

  if (verbosity > 0) {
    std::cout << "Leaving KokkosBlas2_GER()" << std::endl;
  }
}

template <typename tExecSpace>
void run(const blas2_ger_params& params) {
  const auto name = "KokkosBlas2_GER";
  const auto arg_names =
      std::vector<std::string>{"verbosity", "m", "n", "yMode", params.layoutLeft ? "LayoutLeft" : "LayoutRight"};
  const auto args = std::vector<int64_t>{params.verbosity, params.m, params.n, (params.yMode == "transpose"), 1};

  if (params.layoutLeft) {
    if (params.scalarType == "int32") {
      KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas2_GER<std::int32_t, Kokkos::LayoutLeft, tExecSpace>,
                                                 arg_names, args, params.repeat);
    } else if (params.scalarType == "int64") {
      KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas2_GER<std::int64_t, Kokkos::LayoutLeft, tExecSpace>,
                                                 arg_names, args, params.repeat);
    } else if (params.scalarType == "float") {
      KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas2_GER<float, Kokkos::LayoutLeft, tExecSpace>,
                                                 arg_names, args, params.repeat);
    } else if (params.scalarType == "double") {
      KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas2_GER<double, Kokkos::LayoutLeft, tExecSpace>,
                                                 arg_names, args, params.repeat);
    } else if (params.scalarType == "complex_float") {
      KokkosKernelsBenchmark::register_benchmark(
          name, KokkosBlas2_GER<Kokkos::complex<float>, Kokkos::LayoutLeft, tExecSpace>, arg_names, args,
          params.repeat);
    } else {
      KokkosKernelsBenchmark::register_benchmark(
          name, KokkosBlas2_GER<Kokkos::complex<double>, Kokkos::LayoutLeft, tExecSpace>, arg_names, args,
          params.repeat);
    }
  } else {
    if (params.scalarType == "int32") {
      KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas2_GER<std::int32_t, Kokkos::LayoutRight, tExecSpace>,
                                                 arg_names, args, params.repeat);
    } else if (params.scalarType == "int64") {
      KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas2_GER<std::int64_t, Kokkos::LayoutRight, tExecSpace>,
                                                 arg_names, args, params.repeat);
    } else if (params.scalarType == "float") {
      KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas2_GER<float, Kokkos::LayoutRight, tExecSpace>,
                                                 arg_names, args, params.repeat);
    } else if (params.scalarType == "double") {
      KokkosKernelsBenchmark::register_benchmark(name, KokkosBlas2_GER<double, Kokkos::LayoutRight, tExecSpace>,
                                                 arg_names, args, params.repeat);
    } else if (params.scalarType == "complex_float") {
      KokkosKernelsBenchmark::register_benchmark(
          name, KokkosBlas2_GER<Kokkos::complex<float>, Kokkos::LayoutRight, tExecSpace>, arg_names, args,
          params.repeat);
    } else {
      KokkosKernelsBenchmark::register_benchmark(
          name, KokkosBlas2_GER<Kokkos::complex<double>, Kokkos::LayoutRight, tExecSpace>, arg_names, args,
          params.repeat);
    }
  }
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  benchmark::Initialize(&argc, argv);
  benchmark::SetDefaultTimeUnit(benchmark::kSecond);
  KokkosKernelsBenchmark::add_benchmark_context(true);

  const auto params = blas2_ger_params::get_params(argc, argv);

  // std::cout << "In main(): params.repeat = " << params.repeat << std::endl;

  if (params.use_threads) {
#if defined(KOKKOS_ENABLE_THREADS)
    run<Kokkos::Threads>(params);
#else
    std::cout << "ERROR: PThreads requested, but not available.\n";
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
  if (!params.use_cuda && !params.use_hip && !params.use_openmp && !params.use_sycl && !params.use_threads) {
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
