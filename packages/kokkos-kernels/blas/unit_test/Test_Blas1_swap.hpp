#include "KokkosBlas1_swap.hpp"

namespace Test {
namespace Impl {

template <class ScalarType, class DeviceType>
void test_swap(int const vector_length) {
  using execution_space = typename DeviceType::execution_space;
  using memory_space    = typename DeviceType::memory_space;
  using vector_type     = Kokkos::View<ScalarType*, memory_space>;
  using scalar_type     = typename vector_type::non_const_value_type;
  using mag_type        = typename Kokkos::ArithTraits<scalar_type>::mag_type;

  // Note that Xref and Yref need to always be copies of X and Y
  // hence the use of create_mirror instead of create_mirror_view.
  vector_type X("X", vector_length), Y("Y", vector_length);
  typename vector_type::HostMirror Xref = Kokkos::create_mirror(Y);
  typename vector_type::HostMirror Yref = Kokkos::create_mirror(X);

  // Setup values in X, Y and copy them to Xref and Yref
  const scalar_type range = 10 * Kokkos::ArithTraits<scalar_type>::one();
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  Kokkos::fill_random(X, rand_pool, range);
  Kokkos::fill_random(Y, rand_pool, range);

  Kokkos::deep_copy(Xref, Y);
  Kokkos::deep_copy(Yref, X);

  KokkosBlas::swap(X, Y);
  Kokkos::fence();

  typename vector_type::HostMirror Xtest = Kokkos::create_mirror_view(X);
  typename vector_type::HostMirror Ytest = Kokkos::create_mirror_view(Y);
  Kokkos::deep_copy(Xtest, X);
  Kokkos::deep_copy(Ytest, Y);

  const mag_type tol = 10 * Kokkos::ArithTraits<scalar_type>::eps();
  for (int idx = 0; idx < vector_length; ++idx) {
    Test::EXPECT_NEAR_KK_REL(Xtest(idx), Xref(idx), tol);
    Test::EXPECT_NEAR_KK_REL(Ytest(idx), Yref(idx), tol);
  }
}

}  // namespace Impl
}  // namespace Test

template <class ScalarType, class DeviceType>
int test_swap() {
  Test::Impl::test_swap<ScalarType, DeviceType>(0);
  Test::Impl::test_swap<ScalarType, DeviceType>(10);
  Test::Impl::test_swap<ScalarType, DeviceType>(256);
  Test::Impl::test_swap<ScalarType, DeviceType>(1024);

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, swap_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::swap_float");
  test_swap<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, swap_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::swap_double");
  test_swap<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, swap_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::swap_complex_float");
  test_swap<Kokkos::complex<float>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, swap_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::swap_complex_double");
  test_swap<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
