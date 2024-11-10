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

#include "KokkosBlas2_gemv.hpp"
#include <Kokkos_Random.hpp>
#include "KokkosKernels_TestUtils.hpp"

struct Params {
  int use_cuda    = 0;
  int use_hip     = 0;
  int use_openmp  = 0;
  int use_threads = 0;
  int m           = 5000;
  int n           = 5000;
  int repeat      = 1;
  bool layoutLeft = true;
};

void print_options() {
  std::cerr << "Options\n" << std::endl;

  std::cerr << "\tBACKEND: '--threads[numThreads]' | '--openmp [numThreads]' | "
               "'--cuda [cudaDeviceIndex]'"
            << std::endl;
  std::cerr << "\tIf none selected, serial is used." << std::endl;
  std::cerr << "\t[Optional] --repeat :: how many times to repeat overall "
               "spadd (symbolic + repeated numeric)"
            << std::endl;
  std::cerr << "\t[Optional] --layout :: matrix layout ('left' or 'right', "
               "default 'left')"
            << std::endl;
  std::cerr << "\t[Optional] --m      :: number of rows to generate" << std::endl;
  std::cerr << "\t[Optional] --n      :: number of cols to generate" << std::endl;
}

int parse_inputs(Params& params, int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--help") || 0 == Test::string_compare_no_case(argv[i], "-h")) {
      print_options();
      exit(0);  // note: this is before Kokkos::initialize
    } else if (0 == Test::string_compare_no_case(argv[i], "--threads")) {
      params.use_threads = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--openmp")) {
      params.use_openmp = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--cuda")) {
      params.use_cuda = atoi(argv[++i]) + 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--hip")) {
      params.use_hip = atoi(argv[++i]) + 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--layout")) {
      i++;
      if (0 == Test::string_compare_no_case(argv[i], "left"))
        params.layoutLeft = true;
      else if (0 == Test::string_compare_no_case(argv[i], "right"))
        params.layoutLeft = false;
      else {
        std::cerr << "Invalid layout: must be 'left' or 'right'.\n";
        exit(1);
      }
    } else if (0 == Test::string_compare_no_case(argv[i], "--m")) {
      params.m = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--n")) {
      params.n = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--repeat")) {
      // if provided, C will be written to given file.
      // has to have ".bin", or ".crs" extension.
      params.repeat = atoi(argv[++i]);
    } else {
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      print_options();
      return 1;
    }
  }
  return 0;
}

template <typename ExecSpace, typename Layout>
void run(int m, int n, int repeat) {
  // Declare type aliases
  using Scalar   = double;
  using MemSpace = typename ExecSpace::memory_space;
  using Device   = Kokkos::Device<ExecSpace, MemSpace>;

  std::cout << "Running GEMV experiment (" << ExecSpace::name() << ")\n";

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

  // Start timing
  Kokkos::fence();
  Kokkos::Timer timer;
  for (int i = 0; i < repeat; i++) {
    KokkosBlas::gemv("N", 1.0, A, x, 0.0, y);
    ExecSpace().fence();
  }
  // Kokkos Timer set up
  double total = timer.seconds();
  double avg   = total / repeat;
  // Flops calculation
  size_t flopsPerRun = (size_t)2 * m * n;
  printf("Avg GEMV time: %f s.\n", avg);
  printf("Avg GEMV FLOP/s: %.3e\n", flopsPerRun / avg);
}

int main(int argc, char** argv) {
  // Create an instance of Params
  Params params;

  // Argument parsing:
  if (parse_inputs(params, argc, argv)) {
    return 1;
  }
  // const int num_threads = params.use_openmp;
  const int num_threads = std::max(params.use_openmp, params.use_threads);

  const int device_id = std::max(params.use_cuda, params.use_hip) - 1;
  Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(num_threads).set_device_id(device_id));

  // Create booleans to handle pthreads, openmp and cuda params and initialize
  // to true;
  bool useThreads = params.use_threads != 0;
  bool useOMP     = params.use_openmp != 0;
  bool useCUDA    = params.use_cuda != 0;
  bool useHIP     = params.use_hip != 0;

  // Create boolean to handle serial setting if not using open and cuda
  bool useSerial = !useThreads && !useOMP && !useCUDA && !useHIP;

  // Logic for runtime with PThreads
  if (useThreads) {
#if defined(KOKKOS_ENABLE_THREADS)
    if (params.layoutLeft)
      run<Kokkos::Threads, Kokkos::LayoutLeft>(params.m, params.n, params.repeat);
    else
      run<Kokkos::Threads, Kokkos::LayoutRight>(params.m, params.n, params.repeat);
#else
    std::cout << "ERROR:  PThreads requested, but not available.\n";
    return 1;
#endif
  }

  // Logic for runtime with OpenMP
  if (useOMP) {
#if defined(KOKKOS_ENABLE_OPENMP)
    if (params.layoutLeft)
      run<Kokkos::OpenMP, Kokkos::LayoutLeft>(params.m, params.n, params.repeat);
    else
      run<Kokkos::OpenMP, Kokkos::LayoutRight>(params.m, params.n, params.repeat);
#else
    std::cout << "ERROR: OpenMP requested, but not available.\n";
    return 1;
#endif
  }

  // Logic for runtime with Cuda
  if (useCUDA) {
#if defined(KOKKOS_ENABLE_CUDA)
    if (params.layoutLeft)
      run<Kokkos::Cuda, Kokkos::LayoutLeft>(params.m, params.n, params.repeat);
    else
      run<Kokkos::Cuda, Kokkos::LayoutRight>(params.m, params.n, params.repeat);
#else
    std::cout << "ERROR: CUDA requested, but not available.\n";
    return 1;
#endif
  }
  if (useHIP) {
#if defined(KOKKOS_ENABLE_HIP)
    if (params.layoutLeft)
      run<Kokkos::HIP, Kokkos::LayoutLeft>(params.m, params.n, params.repeat);
    else
      run<Kokkos::HIP, Kokkos::LayoutRight>(params.m, params.n, params.repeat);
#else
    std::cout << "ERROR: HIP requested, but not available.\n";
    return 1;
#endif
  }
  // Logic for serial runtime
  if (useSerial) {
#if defined(KOKKOS_ENABLE_SERIAL)
    if (params.layoutLeft)
      run<Kokkos::Serial, Kokkos::LayoutLeft>(params.m, params.n, params.repeat);
    else
      run<Kokkos::Serial, Kokkos::LayoutRight>(params.m, params.n, params.repeat);
#else
    std::cout << "ERROR: Serial device requested, but not available.\n";
    return 1;
#endif
  }
  Kokkos::finalize();
  return 0;
}
