#include <iostream>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_dot.hpp>

int main(int argc, char* argv[]) {
  Kokkos::initialize();
  {
    int N = 100;
    if (argc >= 2) {
      N = atoi(argv[1]);
    }

    Kokkos::View<double*> x("X", N);
    Kokkos::View<double*> y("Y", N);
    Kokkos::deep_copy(x, 3.0);
    Kokkos::deep_copy(y, 2.0);

    double x_y = KokkosBlas::dot(x, y);

    std::cout << "X_dot_Y: " << x_y << " Expected: " << 1.0 * N * (3.0 * 2.0)
              << " Diff: " << x_y - 1.0 * N * (3.0 * 2.0) << std::endl;
  }
  Kokkos::finalize();
}
