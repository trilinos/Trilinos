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

#define KokkosBatched_Test_Trsm_Host_Real 1
#include "KokkosBatched_Test_Trsm_Host.hpp"

using namespace KokkosBatched;

template <typename AlgoTagType>
void run(const int N) {
  typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;

  Kokkos::print_configuration(std::cout, false);

  std::cout << "\n\n Used for Factorization \n\n";

  /// Left, Lower, NoTrans, UnitDiag (used in LU factorization and LU solve)

  PerfTest::Trsm<0, 3, 3, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<0, 5, 5, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<0, 10, 10, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<0, 15, 15, HostSpaceType, AlgoTagType>(N);

  /// Left, Lower, NoTrans, NonUnitDiag

  PerfTest::Trsm<1, 3, 3, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<1, 5, 5, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<1, 10, 10, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<1, 15, 15, HostSpaceType, AlgoTagType>(N);

  /// Right, Upper, NoTrans, UnitDiag

  PerfTest::Trsm<2, 3, 3, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<2, 5, 5, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<2, 10, 10, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<2, 15, 15, HostSpaceType, AlgoTagType>(N);

  /// Right, Upper, NoTrans, NonUnitDiag (used in LU factorization)

  PerfTest::Trsm<3, 3, 3, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<3, 5, 5, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<3, 10, 10, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<3, 15, 15, HostSpaceType, AlgoTagType>(N);

  std::cout << "\n\n Used for Solve \n\n";

  /// Left, Lower, NoTrans, UnitDiag (used in LU solve)

  PerfTest::Trsm<0, 3, 1, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<0, 5, 1, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<0, 10, 1, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<0, 15, 1, HostSpaceType, AlgoTagType>(N);

  /// Left, Upper, Notrans, NonUnitDiag (user in LU solve)

  PerfTest::Trsm<4, 3, 1, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<4, 5, 1, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<4, 10, 1, HostSpaceType, AlgoTagType>(N);
  PerfTest::Trsm<4, 15, 1, HostSpaceType, AlgoTagType>(N);
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

    std::cout << "\n Testing Algo::Trsm::Unblocked\n";
    run<Algo::Trsm::Unblocked>(N);

    std::cout << "\n Testing Algo::Trsm::Blocked\n";
    run<Algo::Trsm::Blocked>(N);

#if defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    std::cout << "\n Testing Algo::Trsm::CompactMKL\n";
    run<Algo::Gemm::CompactMKL>(N);
#endif
  }
#endif
  Kokkos::finalize();

  return 0;
}
