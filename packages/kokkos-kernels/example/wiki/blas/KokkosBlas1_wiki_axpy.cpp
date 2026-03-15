// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <KokkosBlas1_axpby.hpp>

#include <iostream>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    int N = 10;
    if (argc > 1) {
      N = atoi(argv[1]);
    }

    Kokkos::View<double*> x("X", N);
    Kokkos::View<double*> y("Y", N);
    Kokkos::deep_copy(x, 3.0);
    Kokkos::deep_copy(y, 2.0);

    double alpha = 1.5;

    KokkosBlas::axpy(alpha, x, y);

    double sum = 0.0;
    Kokkos::parallel_reduce(
        "CheckValue", N, KOKKOS_LAMBDA(const int& i, double& lsum) { lsum += y(i); }, sum);

    const double expected = N * (2.0 + 1.5 * 3.0);
    std::cout << "Sum: " << sum << ", Expected: " << expected << ", Diff: " << sum - expected << std::endl;
  }

  Kokkos::finalize();
}
