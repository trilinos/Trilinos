#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_scal.hpp>

int main(void) {
  Kokkos::initialize();
  {
    Kokkos::View<double*> x("X", 5);
    auto h_x = Kokkos::create_mirror_view(x);
    Kokkos::deep_copy(h_x, -3.0);
    std::cout << "orignal values:\n"
              << "   {" << h_x(0) << ", " << h_x(1) << ", " << h_x(2) << ", " << h_x(3) << ", " << h_x(4) << "}"
              << std::endl;
    Kokkos::deep_copy(x, h_x);

    KokkosBlas::scal(x, -2.0, x);
    Kokkos::deep_copy(h_x, x);

    std::cout << "final values:\n"
              << "   {" << h_x(0) << ", " << h_x(1) << ", " << h_x(2) << ", " << h_x(3) << ", " << h_x(4) << "}"
              << std::endl;
  }
  Kokkos::finalize();
}
