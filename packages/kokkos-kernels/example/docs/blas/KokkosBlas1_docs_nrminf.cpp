// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_nrminf.hpp>

int main(void) {
  Kokkos::initialize();
  {
    Kokkos::View<double*> x("X", 101);
    Kokkos::parallel_for(
        101, KOKKOS_LAMBDA(const int idx) {
          const double val = static_cast<double>(idx) / 10;
          x(idx)           = (val - 10) * (val - 5) * val + 1.9;
        });

    double x_nrm = KokkosBlas::nrminf(x);

    std::cout << "X_nrm: " << x_nrm << " Expected: " << 50.011 << std::endl;
  }
  Kokkos::finalize();
}
