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

#include <iostream>
#include <chrono>
#include <iomanip>

#include "Kokkos_Random.hpp"
#include "Kokkos_Sort.hpp"

#include "KokkosGraph_Merge.hpp"

struct Result {
  double us;
};

template <typename T, typename Device>
Result bench_merge(size_t aSz, size_t bSz, int nWarmup, int nIters) {
  using execution_space = typename Device::execution_space;

  using view_type = Kokkos::View<T *, Device>;

  view_type a("a", aSz);
  view_type b("b", bSz);
  view_type c("c", a.size() + b.size());

  Kokkos::Random_XorShift64_Pool<execution_space> random_pool(12345);
  fill_random(a, random_pool, 0.0, 100.0);
  fill_random(b, random_pool, 0.0, 100.0);
  Kokkos::sort(a);
  Kokkos::sort(b);
  Kokkos::fence();

  typedef std::chrono::high_resolution_clock Clock;
  typedef std::chrono::duration<double> Duration;

  execution_space space;

  for (int i = 0; i < nWarmup; ++i) {
    KokkosGraph::merge_into(space, c, a, b);
  }
  Kokkos::fence();
  auto start = Clock::now();
  for (int i = 0; i < nIters; ++i) {
    KokkosGraph::merge_into(space, c, a, b);
  }
  Kokkos::fence();
  Duration elapsed(Clock::now() - start);

  Result res;
  res.us = elapsed.count() * 1e6 / nIters;
  return res;
}

int main(int, char **) {
  Kokkos::initialize();

  std::cout << std::setfill(' ');
  std::cout << std::setprecision(3);

#define DO(DEVICE, TYPE, aSz, bSz)                                                            \
  {                                                                                           \
    std::cout << std::setw(9) << aSz << std::setw(9) << bSz << std::setw(9) << #TYPE;         \
    std::cout << std::flush;                                                                  \
    Result res = bench_merge<TYPE, DEVICE>(aSz, bSz, 100, 500);                               \
    std::cout << std::setw(9) << res.us;                                                      \
    std::cout << std::setw(9) << 2 * (aSz + bSz) / (res.us / 1e6);                            \
    std::cout << std::setw(9) << 2 * (aSz + bSz) * sizeof(TYPE) / (res.us / 1e6) / (1 << 30); \
    std::cout << std::endl;                                                                   \
  }

#ifdef KOKKOS_ENABLE_CUDA
  std::cout << "CUDA c=" << Kokkos::Cuda().concurrency() << std::endl;
  std::cout << std::setw(9) << "aSz" << std::setw(9) << "bSz" << std::setw(9) << "type" << std::setw(9) << "us"
            << std::setw(9) << "e/s" << std::setw(9) << "GiB/s" << std::endl;
  DO(Kokkos::Cuda, int32_t, 1e4, 1e4)
  DO(Kokkos::Cuda, int32_t, 1e5, 1e5)
  DO(Kokkos::Cuda, int32_t, 1e6, 1e6)
  DO(Kokkos::Cuda, int32_t, 1e7, 1e7)
  DO(Kokkos::Cuda, int32_t, 2e7, 2e7)
  DO(Kokkos::Cuda, int32_t, 1e4, 2e7)
  DO(Kokkos::Cuda, int32_t, 1e8, 1e8)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
  std::cout << "OpenMP c=" << Kokkos::OpenMP().concurrency() << std::endl;
  std::cout << std::setw(9) << "aSz" << std::setw(9) << "bSz" << std::setw(9) << "type" << std::setw(9) << "us"
            << std::setw(9) << "e/s" << std::setw(9) << "GB/s" << std::endl;
  DO(Kokkos::OpenMP, int32_t, 1e4, 1e4)
  DO(Kokkos::OpenMP, int32_t, 1e5, 1e5)
  DO(Kokkos::OpenMP, int32_t, 1e6, 1e6)
  DO(Kokkos::OpenMP, int32_t, 1e7, 1e7)
  DO(Kokkos::OpenMP, int32_t, 2e7, 2e7)
  DO(Kokkos::OpenMP, int32_t, 1e4, 4e7)
#endif

  Kokkos::finalize();

  return 0;
}
