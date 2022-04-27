/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>
#include <blas/KokkosBlas1_dot.hpp>
#include <Kokkos_Random.hpp>

// For RPS implementation
#include "KokkosBlas_dot_perf_test.hpp"
#include "KokkosKernels_TestUtils.hpp"

struct Params {
  int use_cuda    = 0;
  int use_openmp  = 0;
  int use_threads = 0;
  // m is vector length
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
  std::cerr << "\t[Optional] --m      :: desired length of test vectors; test "
               "vectors will have the same length"
            << std::endl;
}

int parse_inputs(Params& params, int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--help") ||
        0 == Test::string_compare_no_case(argv[i], "-h")) {
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
      std::cerr << "Unrecognized command line argument #" << i << ": "
                << argv[i] << std::endl;
      print_options();
      return 1;
    }
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// The Level 1 BLAS perform scalar, vector and vector-vector operations;
//
// https://github.com/kokkos/kokkos-kernels/wiki/BLAS-1%3A%3Adot
//
// Usage: result = KokkosBlas::dot(x,y); KokkosBlas::dot(r,x,y);
// Multiplies each value of x(i) [x(i,j)] with y(i) or [y(i,j)] and computes the
// sum. (If x and y have scalar type Kokkos::complex, the complex conjugate of
// x(i) or x(i,j) will be used.) VectorX: A rank-1 Kokkos::View VectorY: A
// rank-1 Kokkos::View ReturnVector: A rank-0 or rank-1 Kokkos::View
//
// REQUIREMENTS:
// Y.rank == 1 or X.rank == 1
// Y.extent(0) == X.extent(0)

// Dot Test design:
// 1) create 1D View containing 1D matrix, aka a vector; this will be your X
// input matrix; 2) create 1D View containing 1D matrix, aka a vector; this will
// be your Y input matrix; 3) perform the dot operation on the two inputs, and
// capture result in "result"

// Here, m represents the desired length for each 1D matrix;
// "m" is used here, because code from another test was adapted for this test.
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class ExecSpace>
void run(int m, int repeat) {
  // Declare type aliases
  using Scalar   = double;
  using MemSpace = typename ExecSpace::memory_space;
  using Device   = Kokkos::Device<ExecSpace, MemSpace>;

  std::cout << "Running BLAS Level 1 DOT perfomrance experiment ("
            << ExecSpace::name() << ")\n";

  std::cout << "Each test input vector has a length of " << m << std::endl;

  // Create 1D view w/ Device as the ExecSpace; this is an input vector
  // A(view_alloc(WithoutInitializing, "label"), m, n);
  Kokkos::View<Scalar*, Device> x(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "x"), m);

  // Create 1D view w/ Device as the ExecSpace; this is the output vector
  Kokkos::View<Scalar*, Device> y(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "y"), m);

  // Declaring variable pool w/ a seeded random number;
  // a parallel random number generator, so you
  // won't get the same number with a given seed each time
  Kokkos::Random_XorShift64_Pool<ExecSpace> pool(123);

  Kokkos::fill_random(x, pool, 10.0);
  Kokkos::fill_random(y, pool, 10.0);

  // do a warm up run of dot:
  KokkosBlas::dot(x, y);

  // The live test of dot:

  Kokkos::fence();
  Kokkos::Timer timer;

  for (int i = 0; i < repeat; i++) {
    KokkosBlas::dot(x, y);
    ExecSpace().fence();
  }

  // Kokkos Timer set up
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

  Kokkos::initialize(Kokkos::InitArguments(num_threads, -1, device_id));

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
    std::cout << "ERROR: Serial device requested, but not available.\n";
    return 1;
#endif
  }
  Kokkos::finalize();
  return 0;
}
