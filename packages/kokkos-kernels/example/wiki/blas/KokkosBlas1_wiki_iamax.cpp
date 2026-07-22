#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas1_iamax.hpp>

int main(int argc, char* argv[]) {
  Kokkos::initialize();
  {
    int N = 100;
    if (argc >= 2) {
      N = atoi(argv[1]);
    }

    using ViewType  = Kokkos::View<double*>;
    using Scalar    = typename ViewType::non_const_value_type;
    using AT        = KokkosKernels::ArithTraits<Scalar>;
    using mag_type  = typename AT::mag_type;
    using size_type = typename ViewType::size_type;

    ViewType x("X", N);

    typename ViewType::host_mirror_type h_x = Kokkos::create_mirror_view(x);

    Kokkos::Random_XorShift64_Pool<typename ViewType::device_type::execution_space> rand_pool(13718);
    Kokkos::fill_random(x, rand_pool, Scalar(10));

    Kokkos::deep_copy(h_x, x);

    size_type max_loc = KokkosBlas::iamax(x);

    mag_type expected_result   = KokkosKernels::ArithTraits<mag_type>::min();
    size_type expected_max_loc = 0;
    for (int i = 0; i < N; i++) {
      mag_type val = AT::abs(h_x(i));
      if (val > expected_result) {
        expected_result  = val;
        expected_max_loc = i + 1;
      }
    }

    std::cout << "Iamax of X: " << max_loc << ", Expected: " << expected_max_loc << std::endl;
  }
  Kokkos::finalize();
}
