// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_fill.hpp>

int main(void) {
  Kokkos::initialize();
  {
    Kokkos::View<double*> x("x", 5);
    auto h_x = Kokkos::create_mirror_view(x);
    h_x(0)   = 0.0;
    h_x(1)   = -1.0;
    h_x(2)   = -2.0;
    h_x(3)   = 3.0;
    h_x(4)   = -4.0;

    std::cout << "orignal values:\n"
              << "   x: {" << h_x(0) << ", " << h_x(1) << ", " << h_x(2) << ", " << h_x(3) << ", " << h_x(4) << "}"
              << std::endl;
    Kokkos::deep_copy(x, h_x);

    KokkosBlas::fill(x, 7);

    Kokkos::deep_copy(h_x, x);
    std::cout << "final values:\n"
              << "   x: {" << h_x(0) << ", " << h_x(1) << ", " << h_x(2) << ", " << h_x(3) << ", " << h_x(4) << "}"
              << std::endl;
  }
  Kokkos::finalize();
}
