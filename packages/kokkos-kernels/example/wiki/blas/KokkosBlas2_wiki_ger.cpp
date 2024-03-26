#include <Kokkos_Core.hpp>
#include <KokkosBlas2_ger.hpp>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    constexpr int M = 5;
    constexpr int N = 4;

    Kokkos::View<double**> A("A", M, N);
    Kokkos::View<double*> x("X", M);
    Kokkos::View<double*> y("Y", N);

    Kokkos::deep_copy(A, 1.0);
    Kokkos::deep_copy(x, 3.0);
    Kokkos::deep_copy(y, 1.3);

    const double alpha = Kokkos::ArithTraits<double>::one();

    KokkosBlas::ger("T", alpha, x, y, A);
  }
  Kokkos::finalize();
}
