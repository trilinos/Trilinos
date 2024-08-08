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
#include <KokkosBlas1_team_dot.hpp>
#include <Kokkos_Random.hpp>
#include "KokkosKernels_TestUtils.hpp"

struct Params {
  int use_cuda    = 0;
  int use_openmp  = 0;
  int use_threads = 0;
  // m is vector length, or number of rows
  int m      = 100000;
  int repeat = 1;
};

void print_options() {
  std::cerr << "Options:\n" << std::endl;

  std::cerr << "\tBACKEND: '--threads[numThreads]' | '--openmp [numThreads]' | "
               "'--cuda [cudaDeviceIndex]'"
            << std::endl;
  std::cerr << "\tIf no BACKEND selected, serial is the default." << std::endl;
  std::cerr << "\t[Optional] --repeat :: how many times to repeat overall "
               "dot (symbolic + repeated numeric)"
            << std::endl;
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
    } else if (0 == Test::string_compare_no_case(argv[i], "--m")) {
      params.m = atoi(argv[++i]);
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
void run(int m, int repeat) {
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

  std::cout << "Running BLAS Level 1 Kokkos Teams-based implementation DOT "
               "performance experiment ("
            << ExecSpace::name() << ")\n";

  std::cout << "Each test input vector has a length of " << m << std::endl;

  // Warm up run of dot:
  teamDotFunctor<Kokkos::View<Scalar*, MemSpace>, ExecSpace> teamDotFunctorWarmUpInstance(x, y);

  Kokkos::parallel_for("TeamDotUsage -- Warm Up Run", policy(1, Kokkos::AUTO), teamDotFunctorWarmUpInstance);

  Kokkos::fence();
  Kokkos::Timer timer;

  // Live test of dot:

  teamDotFunctor<Kokkos::View<Scalar*, MemSpace>, ExecSpace> teamDotFunctorLiveTestInstance(x, y);
  Kokkos::parallel_for("TeamDotUsage -- Live Test", policy(1, Kokkos::AUTO), teamDotFunctorLiveTestInstance);

  ExecSpace().fence();

  // Kokkos Timer set up and data capture
  double total = timer.seconds();
  double avg   = total / repeat;
  // Flops calculation for a 1D matrix dot product per test run;
  size_t flopsPerRun = (size_t)2 * m;
  printf("Avg DOT time: %f s.\n", avg);
  printf("Avg DOT FLOP/s: %.3e\n", flopsPerRun / avg);
}

int main(int argc, char** argv) {
  Params params;

  if (parse_inputs(params, argc, argv)) {
    return 1;
  }

  const int device_id = params.use_cuda - 1;

  const int num_threads = std::max(params.use_openmp, params.use_threads);

  Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(num_threads).set_device_id(device_id));

  bool useThreads = params.use_threads != 0;
  bool useOMP     = params.use_openmp != 0;
  bool useCUDA    = params.use_cuda != 0;
  bool useSerial  = !useThreads && !useOMP && !useCUDA;

  if (useThreads) {
#if defined(KOKKOS_ENABLE_THREADS)
    run<Kokkos::Threads>(params.m, params.repeat);
#else
    std::cout << "ERROR:  PThreads requested, but not available.\n";
    return 1;
#endif
  }

  if (useOMP) {
#if defined(KOKKOS_ENABLE_OPENMP)
    run<Kokkos::OpenMP>(params.m, params.repeat);
#else
    std::cout << "ERROR: OpenMP requested, but not available.\n";
    return 1;
#endif
  }

  if (useCUDA) {
#if defined(KOKKOS_ENABLE_CUDA)
    run<Kokkos::Cuda>(params.m, params.repeat);
#else
    std::cout << "ERROR: CUDA requested, but not available.\n";
    return 1;
#endif
  }
  if (useSerial) {
#if defined(KOKKOS_ENABLE_SERIAL)
    run<Kokkos::Serial>(params.m, params.repeat);
#else
    std::cout << "ERROR: Serial device requested, but not available; here, "
                 "implementation of dot is explicitly parallel.\n";
    return 1;
#endif
  }
  Kokkos::finalize();
  return 0;
}
