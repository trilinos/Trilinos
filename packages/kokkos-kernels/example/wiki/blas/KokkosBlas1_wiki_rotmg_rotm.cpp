#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include "KokkosBlas1_rotmg.hpp"
#include "KokkosBlas1_rotm.hpp"
#include "KokkosKernels_PrintUtils.hpp"

int main(int argc, char* argv[]) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using Scalar          = double;
  using Vector          = Kokkos::View<Scalar*, execution_space>;
  using ParamView       = Kokkos::View<Scalar[5], execution_space>;
  using ScalarView      = Kokkos::View<Scalar, execution_space>;

  Kokkos::initialize(argc, argv);
  {
    const int N = 10;
    Vector x("x", N);
    Vector y("y", N);
    ScalarView d1("d1");
    ScalarView d2("d2");
    ParamView param("param");

    // Populate x,y with uniform random values between 0 and 10
    Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
    Kokkos::fill_random(x, rand_pool, Scalar(10));
    Kokkos::fill_random(y, rand_pool, Scalar(10));

    // Populate input vector scaling factors with 1
    Kokkos::deep_copy(d1, Scalar(1));
    Kokkos::deep_copy(d2, Scalar(1));

    std::cout << "x,y before applying modified Givens rotation:\n";
    KokkosKernels::Impl::kk_print_1Dview(std::cout, x);
    KokkosKernels::Impl::kk_print_1Dview(std::cout, y);

    // Calculate Givens rotation coefficients to eliminate y(0)
    KokkosBlas::rotmg<execution_space, ScalarView, ScalarView, ParamView>(
        execution_space(), d1, d2, Kokkos::subview(x, 0), Kokkos::subview(y, 0), param);

    auto paramHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), param);

    std::cout << "\nrotmg output (rotation parameters) to eliminate y(0):\n";
    std::cout << "d1 = ";
    KokkosKernels::Impl::kk_print_1Dview(std::cout, d1);
    std::cout << "d2 = ";
    KokkosKernels::Impl::kk_print_1Dview(std::cout, d2);
    std::cout << "flag = " << paramHost(0) << '\n';
    std::cout << "h components = ";
    for (int i = 0; i < 4; i++) std::cout << paramHost(1 + i) << " ";
    std::cout << '\n';

    // Zero out y(0), which was left unmodified by rotmg.
    Kokkos::deep_copy(Kokkos::subview(y, 0), Scalar(0));

    // Apply the rotation to the remaining entries of x and y
    KokkosBlas::rotm(execution_space(), Kokkos::subview(x, Kokkos::make_pair(1, N)),
                     Kokkos::subview(y, Kokkos::make_pair(1, N)), param);

    // Apply scaling factors: sqrt(d1) and sqrt(d2) to x and y respectively
    Kokkos::parallel_for(
        Kokkos::RangePolicy<execution_space>(0, N), KOKKOS_LAMBDA(int i) {
          x(i) *= Kokkos::sqrt(d1());
          y(i) *= Kokkos::sqrt(d2());
        });

    std::cout << "\nx,y after applying modified Givens rotation and scaling by [sqrt(d1), sqrt(d2)]):\n";
    KokkosKernels::Impl::kk_print_1Dview(std::cout, x);
    KokkosKernels::Impl::kk_print_1Dview(std::cout, y);
  }
  Kokkos::finalize();
}
