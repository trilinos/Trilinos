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

#include "KokkosBlas3_gemm.hpp"
#include <Kokkos_Random.hpp>
#include "KokkosKernels_TestUtils.hpp"

struct Params {
  int use_cuda    = 0;
  int use_openmp  = 0;
  int use_threads = 0;
  int m           = 1000;
  int n           = 1000;
  int k           = 1000;
  int repeat      = 1;
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
  std::cerr << "\t[Optional] --m      :: Rows in A" << std::endl;
  std::cerr << "\t[Optional] --n      :: Columns in A / Rows in B" << std::endl;
  std::cerr << "\t[Optional] --k      :: Columns in B" << std::endl;
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
    } else if (0 == Test::string_compare_no_case(argv[i], "--n")) {
      params.n = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--k")) {
      params.k = atoi(argv[++i]);
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

template <typename ExecSpace, typename ALayout, typename BLayout>
void runImpl(int m, int n, int k, int repeat) {
  using Scalar   = double;
  using MemSpace = typename ExecSpace::memory_space;
  using Device   = Kokkos::Device<ExecSpace, MemSpace>;
  Kokkos::View<Scalar**, ALayout, Device> A(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "A"), m, n);
  Kokkos::View<Scalar**, BLayout, Device> B(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "B"), n, k);
  Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device> C(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "C"), m, k);
  Kokkos::Random_XorShift64_Pool<ExecSpace> pool(123);
  Kokkos::fill_random(A, pool, 10.0);
  Kokkos::fill_random(B, pool, 10.0);
  // Do a warm-up run
  KokkosBlas::gemm("N", "N", 1.0, A, B, 0.0, C);
  // Now, start timing
  Kokkos::fence();
  Kokkos::Timer timer;
  for (int i = 0; i < repeat; i++) {
    KokkosBlas::gemm("N", "N", 1.0, A, B, 0.0, C);
    ExecSpace().fence();
  }
  double total       = timer.seconds();
  double avg         = total / repeat;
  size_t flopsPerRun = (size_t)2 * m * n * k;
  printf("Avg GEMM FLOP/s: %.3e --- Avg time: %f\n", flopsPerRun / avg, avg);
}

template <typename ExecSpace>
void run(int m, int n, int k, int repeat) {
  using LL = Kokkos::LayoutLeft;
  using LR = Kokkos::LayoutRight;
  std::cout << "** Running GEMM experiments (" << ExecSpace::name() << ") **\n";
  std::cout << "Running: A LayoutLeft, B LayoutLeft  : ";
  runImpl<ExecSpace, LL, LL>(m, n, k, repeat);
  std::cout << "Running: A LayoutLeft, B LayoutRight : ";
  runImpl<ExecSpace, LL, LR>(m, n, k, repeat);
  std::cout << "Running: A LayoutRight, B LayoutLeft : ";
  runImpl<ExecSpace, LR, LL>(m, n, k, repeat);
  std::cout << "Running: A LayoutRight, B LayoutRight: ";
  runImpl<ExecSpace, LR, LR>(m, n, k, repeat);
}

int main(int argc, char** argv) {
  Params params;

  if (parse_inputs(params, argc, argv)) {
    return 1;
  }
  const int num_threads =
      params.use_openmp;  // Assumption is that use_openmp variable is provided
                          // as number of threads
  const int device_id = params.use_cuda - 1;

  Kokkos::initialize(Kokkos::InitArguments(num_threads, -1, device_id));

  bool useOMP  = params.use_openmp != 0;
  bool useCUDA = params.use_cuda != 0;

  bool useSerial = !useOMP && !useCUDA;

  if (useOMP) {
#if defined(KOKKOS_ENABLE_OPENMP)
    run<Kokkos::OpenMP>(params.m, params.n, params.k, params.repeat);
#else
    std::cout << "ERROR: OpenMP requested, but not available.\n";
    return 1;
#endif
  }
  if (useCUDA) {
#if defined(KOKKOS_ENABLE_CUDA)
    run<Kokkos::Cuda>(params.m, params.n, params.k, params.repeat);
#else
    std::cout << "ERROR: CUDA requested, but not available.\n";
    return 1;
#endif
  }
  if (useSerial) {
#if defined(KOKKOS_ENABLE_SERIAL)
    run<Kokkos::Serial>(params.m, params.n, params.k, params.repeat);
#else
    std::cout << "ERROR: Serial device requested, but not available.\n";
    return 1;
#endif
  }
  Kokkos::finalize();
  return 0;
}
