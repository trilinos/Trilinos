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
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#define KokkosBatched_Test_LU_Host_Real 1
#include "KokkosBatched_Test_LU_Host.hpp"

using namespace KokkosBatched;

template <typename AlgoTagType>
void run(const int N) {
  typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;

  Kokkos::print_configuration(std::cout, false);

  PerfTest::LU<3, HostSpaceType, AlgoTagType>(N);
  PerfTest::LU<5, HostSpaceType, AlgoTagType>(N);
  PerfTest::LU<10, HostSpaceType, AlgoTagType>(N);
  PerfTest::LU<15, HostSpaceType, AlgoTagType>(N);
}

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);

#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
  int N = 128 * 128;

  for (int i = 1; i < argc; ++i) {
    const std::string& token = argv[i];
    if (token == std::string("-N")) N = std::atoi(argv[++i]);
  }

  {
    std::cout << " N = " << N << std::endl;

    std::cout << "\n Testing Algo::LU::Unblocked\n";
    run<Algo::LU::Unblocked>(N);

    std::cout << "\n Testing Algo::LU::Blocked\n";
    run<Algo::LU::Blocked>(N);

#if defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    std::cout << "\n Testing Algo::LU::CompactMKL\n";
    run<Algo::Gemm::CompactMKL>(N);
#endif
  }
#endif

  Kokkos::finalize();

  return 0;
}
