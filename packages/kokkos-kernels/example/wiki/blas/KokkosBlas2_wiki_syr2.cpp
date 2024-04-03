#include <Kokkos_Core.hpp>
#include <KokkosBlas2_syr2.hpp>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    constexpr int M = 5;

    Kokkos::View<double**> A("A", M, M);
    Kokkos::View<double*> x("X", M);
    Kokkos::View<double*> y("Y", M);

    Kokkos::deep_copy(A, 1.0);
    Kokkos::deep_copy(x, 3.0);
    Kokkos::deep_copy(y, 1.3);

    const double alpha = double(1.0);

    KokkosBlas::syr2("T", "U", alpha, x, y, A);
  }
  Kokkos::finalize();
}
