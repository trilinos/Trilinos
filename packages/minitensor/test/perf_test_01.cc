// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <chrono>
#include <vector>

#include "MiniTensor.h"

struct Now {
  typedef std::chrono::time_point<std::chrono::high_resolution_clock> Impl;
  Impl impl;
};

Now now() {
  Now t;
  t.impl = std::chrono::high_resolution_clock::now();
  return t;
}

double operator-(Now b, Now a) {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(b.impl - a.impl)
             .count() *
         1e-9;
}

int
main(int ac, char* av[])
{
  Kokkos::initialize();

  int const n = 1000 * 1000;
  Kokkos::View<double*, Kokkos::DefaultExecutionSpace> results("results", n);
  auto t0 = now();
  for (int repeat_i = 0; repeat_i < 20; ++repeat_i) {
    Kokkos::parallel_for(Kokkos::RangePolicy<int, Kokkos::DefaultHostExecutionSpace>(0, n), [=](int i) {
      minitensor::Tensor<double, 3> B;
      B(0,0) = 3.14 * i;
      B(1,1) = 3.14 - double(i);
      B(2,2) = 3.14 + double(i);
      results(i) = trace(B);
    });
  }
  auto t1 = now();

  std::cout << "result(4242) " << results(4242) << '\n';
  std::cout << "total time " << (t1 - t0) << '\n';

  Kokkos::finalize();
}
