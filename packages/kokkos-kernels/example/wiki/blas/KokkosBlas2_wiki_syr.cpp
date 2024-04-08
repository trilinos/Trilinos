#include <Kokkos_Core.hpp>
#include <KokkosBlas2_syr.hpp>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    constexpr int M = 5;

    Kokkos::View<double**> A("A", M, M);
    Kokkos::View<double*> x("X", M);

    Kokkos::deep_copy(A, 1.0);
    Kokkos::deep_copy(x, 3.0);

    const double alpha = double(1.0);

    KokkosBlas::syr("T", "U", alpha, x, A);
  }
  Kokkos::finalize();
}
