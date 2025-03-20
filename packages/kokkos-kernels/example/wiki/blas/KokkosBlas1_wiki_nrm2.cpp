#include <cmath>
#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_nrm2.hpp>

int main(void) {
  Kokkos::initialize();
  {
    Kokkos::View<double*> x("X", 100);
    Kokkos::deep_copy(x, 3.0);

    double x_nrm = KokkosBlas::nrm2(x);

    std::cout << "X_nrm: " << x_nrm << " Expected: " << std::sqrt(100 * 3.0 * 3.0) << std::endl;
  }
  Kokkos::finalize();
}
