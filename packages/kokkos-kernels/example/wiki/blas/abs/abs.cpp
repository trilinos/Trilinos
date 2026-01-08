// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
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
