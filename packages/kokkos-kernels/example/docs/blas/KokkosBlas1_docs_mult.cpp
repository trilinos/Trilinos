// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_mult.hpp>

int main(void) {
  Kokkos::initialize();
  {
    Kokkos::View<double*> a("a", 5), x("x", 5), y("y", 5);
    auto h_y = Kokkos::create_mirror_view(y);
    h_y(0)   = 3.0;
    h_y(1)   = 2.0;
    h_y(2)   = 3.0;
    h_y(3)   = -2.0;
    h_y(4)   = 4.0;
    auto h_a = Kokkos::create_mirror_view(a);
    h_a(0)   = -3.0;
    h_a(1)   = -1.0;
    h_a(2)   = 2.0;
    h_a(3)   = -2.0;
    h_a(4)   = -3.0;
    auto h_x = Kokkos::create_mirror_view(x);
    h_x(0)   = 0.0;
    h_x(1)   = -1.0;
    h_x(2)   = -2.0;
    h_x(3)   = 3.0;
    h_x(4)   = -4.0;
    Kokkos::deep_copy(y, h_y);
    Kokkos::deep_copy(a, h_a);
    Kokkos::deep_copy(x, h_x);

    std::cout << "orignal values:\n"
              << "   beta: 2.0,  alpha: -1.0\n"
              << "   y: {" << h_y(0) << ", " << h_y(1) << ", " << h_y(2) << ", " << h_y(3) << ", " << h_y(4) << "}\n"
              << "   a: {" << h_a(0) << ", " << h_a(1) << ", " << h_a(2) << ", " << h_a(3) << ", " << h_a(4) << "}\n"
              << "   x: {" << h_x(0) << ", " << h_x(1) << ", " << h_x(2) << ", " << h_x(3) << ", " << h_x(4) << "}\n"
              << std::endl;

    KokkosBlas::mult(2.0, y, -1.0, a, x);

    Kokkos::deep_copy(h_y, y);
    std::cout << "final values:\n"
              << "   y: {" << h_y(0) << ", " << h_y(1) << ", " << h_y(2) << ", " << h_y(3) << ", " << h_y(4) << "}"
              << std::endl;
  }
  Kokkos::finalize();
}
