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
#include <KokkosBlas1_abs.hpp>

int main(int argc, char* argv[]) {
  Kokkos::initialize();

  int N = atoi(argv[1]);

  Kokkos::View<double*> x("X", N);
  Kokkos::View<double*> y("Y", N);
  Kokkos::deep_copy(x, -1.0);

  KokkosBlas::abs(y, x);

  double sum = 0.0;
  Kokkos::parallel_reduce(
      "CheckValue", N, KOKKOS_LAMBDA(const int& i, double& lsum) { lsum += y(i); }, sum);

  printf("Sum: %lf Expected: %lf Diff: %e\n", sum, 1.0 * N, sum - 1.0 * N);

  Kokkos::finalize();
}
