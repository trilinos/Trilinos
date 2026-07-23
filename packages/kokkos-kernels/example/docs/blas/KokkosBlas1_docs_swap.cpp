#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_swap.hpp>

int main(void) {
  Kokkos::initialize();
  {
    Kokkos::View<double*> x("x", 5), y("y", 5);
    auto h_x = Kokkos::create_mirror_view(x);
    auto h_y = Kokkos::create_mirror_view(y);
    h_x(0)   = 0.0;
    h_x(1)   = 1.0;
    h_x(2)   = 2.0;
    h_x(3)   = 3.0;
    h_x(4)   = 4.0;
    h_y(0)   = 9.0;
    h_y(1)   = 8.0;
    h_y(2)   = 7.0;
    h_y(3)   = 6.0;
    h_y(4)   = 5.0;

    std::cout << "orignal values:\n"
              << "   x: {" << h_x(0) << ", " << h_x(1) << ", " << h_x(2) << ", " << h_x(3) << ", " << h_x(4) << "}\n"
              << "   y: {" << h_y(0) << ", " << h_y(1) << ", " << h_y(2) << ", " << h_y(3) << ", " << h_y(4) << "}"
              << std::endl;
    Kokkos::deep_copy(x, h_x);
    Kokkos::deep_copy(y, h_y);

    KokkosBlas::swap(x, y);
    Kokkos::deep_copy(h_x, x);
    Kokkos::deep_copy(h_y, y);

    std::cout << "final values:\n"
              << "   x: {" << h_x(0) << ", " << h_x(1) << ", " << h_x(2) << ", " << h_x(3) << ", " << h_x(4) << "}\n"
              << "   y: {" << h_y(0) << ", " << h_y(1) << ", " << h_y(2) << ", " << h_y(3) << ", " << h_y(4) << "}"
              << std::endl;
  }
  Kokkos::finalize();
}
