// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_update.hpp>

int main(void) {
  Kokkos::initialize();
  {
    Kokkos::View<double*> x("x", 5), y("y", 5), z("z", 5);
    auto h_x = Kokkos::create_mirror_view(x);
    h_x(0)   = 0.0;
    h_x(1)   = -1.0;
    h_x(2)   = -2.0;
    h_x(3)   = 3.0;
    h_x(4)   = -4.0;
    auto h_y = Kokkos::create_mirror_view(y);
    h_y(0)   = 3.0;
    h_y(1)   = 2.0;
    h_y(2)   = 3.0;
    h_y(3)   = -2.0;
    h_y(4)   = 4.0;
    auto h_z = Kokkos::create_mirror_view(z);
    h_z(0)   = -3.0;
    h_z(1)   = -1.0;
    h_z(2)   = 2.0;
    h_z(3)   = -2.0;
    h_z(4)   = -3.0;
    Kokkos::deep_copy(x, h_x);
    Kokkos::deep_copy(y, h_y);
    Kokkos::deep_copy(z, h_z);

    std::cout << "orignal values:\n"
              << "   alpha: -1.0,  beta: 2.0,  gamma: \n"
              << "   x: {" << h_x(0) << ", " << h_x(1) << ", " << h_x(2) << ", " << h_x(3) << ", " << h_x(4) << "}\n"
              << "   y: {" << h_y(0) << ", " << h_y(1) << ", " << h_y(2) << ", " << h_y(3) << ", " << h_y(4) << "}\n"
              << "   z: {" << h_z(0) << ", " << h_z(1) << ", " << h_z(2) << ", " << h_z(3) << ", " << h_z(4) << "}\n"
              << std::endl;

    KokkosBlas::update(-1.0, x, 2.0, y, 1.0, z);

    Kokkos::deep_copy(h_z, z);
    std::cout << "final values:\n"
              << "   z: {" << h_z(0) << ", " << h_z(1) << ", " << h_z(2) << ", " << h_z(3) << ", " << h_z(4) << "}"
              << std::endl;
  }
  Kokkos::finalize();
}
