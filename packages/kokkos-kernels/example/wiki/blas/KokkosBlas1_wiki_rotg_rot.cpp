#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include "KokkosBlas1_rotg.hpp"
#include "KokkosBlas1_rot.hpp"
#include "KokkosKernels_PrintUtils.hpp"

int main(int argc, char* argv[]) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using Scalar          = double;
  using Vector          = Kokkos::View<Scalar*, execution_space>;
  using ScalarView      = Kokkos::View<Scalar, execution_space>;

  Kokkos::initialize(argc, argv);
  {
    const int N = 10;
    Vector x("x", N);
    Vector y("y", N);

    // Populate x,y with uniform random values between 0 and 10
    Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
    Kokkos::fill_random(x, rand_pool, Scalar(10));
    Kokkos::fill_random(y, rand_pool, Scalar(10));

    std::cout << "x,y before applying Givens rotation:\n";
    KokkosKernels::Impl::kk_print_1Dview(std::cout, x);
    KokkosKernels::Impl::kk_print_1Dview(std::cout, y);

    ScalarView c("c");
    ScalarView s("s");

    // Calculate Givens rotation coefficients to eliminate y(0)
    KokkosBlas::rotg<execution_space, ScalarView, ScalarView>(execution_space(), Kokkos::subview(x, 0),
                                                              Kokkos::subview(y, 0), c, s);

    std::cout << "\nrotg output (rotation parameters) to eliminate y(0):\n";
    std::cout << "c = ";
    KokkosKernels::Impl::kk_print_1Dview(std::cout, c);
    std::cout << "s = ";
    KokkosKernels::Impl::kk_print_1Dview(std::cout, s);
    std::cout << "r = x(0) = ";
    KokkosKernels::Impl::kk_print_1Dview(std::cout, Kokkos::subview(x, 0));
    std::cout << "z = ";
    KokkosKernels::Impl::kk_print_1Dview(std::cout, Kokkos::subview(y, 0));

    // Zero out y(0), which now contains the output parameter z.
    // This completes the replacement of [x(0), y(0)] with [r, 0].
    Kokkos::deep_copy(Kokkos::subview(y, 0), Scalar(0));

    // Apply the rotation to the remaining entries of x and y
    KokkosBlas::rot(execution_space(), Kokkos::subview(x, Kokkos::make_pair(1, N)),
                    Kokkos::subview(y, Kokkos::make_pair(1, N)), c, s);

    std::cout << "\nx,y after applying Givens rotation:\n";
    KokkosKernels::Impl::kk_print_1Dview(std::cout, x);
    KokkosKernels::Impl::kk_print_1Dview(std::cout, y);
  }
  Kokkos::finalize();
}
