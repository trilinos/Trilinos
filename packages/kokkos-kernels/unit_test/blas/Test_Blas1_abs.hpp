#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas1_abs.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_abs(int N) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef Kokkos::Details::ArithTraits<ScalarA> AT;

  typedef Kokkos::View<
      ScalarA * [2],
      typename std::conditional<std::is_same<typename ViewTypeA::array_layout,
                                             Kokkos::LayoutStride>::value,
                                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,
      Device>
      BaseTypeA;
  typedef Kokkos::View<
      ScalarB * [2],
      typename std::conditional<std::is_same<typename ViewTypeB::array_layout,
                                             Kokkos::LayoutStride>::value,
                                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,
      Device>
      BaseTypeB;

  typename AT::mag_type eps = AT::epsilon() * 10;

  BaseTypeA b_x("X", N);
  BaseTypeB b_y("Y", N);
  BaseTypeB b_org_y("Org_Y", N);

  ViewTypeA x                        = Kokkos::subview(b_x, Kokkos::ALL(), 0);
  ViewTypeB y                        = Kokkos::subview(b_y, Kokkos::ALL(), 0);
  typename ViewTypeA::const_type c_x = x;

  typename BaseTypeA::HostMirror h_b_x = Kokkos::create_mirror_view(b_x);
  typename BaseTypeB::HostMirror h_b_y = Kokkos::create_mirror_view(b_y);

  typename ViewTypeA::HostMirror h_x = Kokkos::subview(h_b_x, Kokkos::ALL(), 0);
  typename ViewTypeB::HostMirror h_y = Kokkos::subview(h_b_y, Kokkos::ALL(), 0);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(b_x, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(b_y, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(b_org_y, b_y);

  Kokkos::deep_copy(h_b_x, b_x);
  Kokkos::deep_copy(h_b_y, b_y);

  // Run with nonconst input
  KokkosBlas::abs(y, x);
  // Copy result to host (h_y is subview of h_b_y)
  Kokkos::deep_copy(h_b_y, b_y);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(h_y(i), AT::abs(h_x(i)), eps * AT::abs(h_x(i)));
  }
  // Run with const input
  // Reset output
  Kokkos::deep_copy(b_y, b_org_y);
  KokkosBlas::abs(y, c_x);
  Kokkos::deep_copy(h_b_y, b_y);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(h_y(i), AT::abs(h_x(i)), eps * AT::abs(h_x(i)));
  }
}

template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_abs_mv(int N, int K) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef Kokkos::Details::ArithTraits<ScalarA> AT;

  typedef multivector_layout_adapter<ViewTypeA> vfA_type;
  typedef multivector_layout_adapter<ViewTypeB> vfB_type;

  typename vfA_type::BaseType b_x("A", N, K);
  typename vfB_type::BaseType b_y("B", N, K);
  typename vfB_type::BaseType b_org_y("B", N, K);

  ViewTypeA x = vfA_type::view(b_x);
  ViewTypeB y = vfB_type::view(b_y);

  typedef multivector_layout_adapter<typename ViewTypeA::HostMirror> h_vfA_type;
  typedef multivector_layout_adapter<typename ViewTypeB::HostMirror> h_vfB_type;

  typename h_vfA_type::BaseType h_b_x = Kokkos::create_mirror_view(b_x);
  typename h_vfB_type::BaseType h_b_y = Kokkos::create_mirror_view(b_y);

  typename ViewTypeA::HostMirror h_x = h_vfA_type::view(h_b_x);
  typename ViewTypeB::HostMirror h_y = h_vfB_type::view(h_b_y);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(b_x, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(b_y, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(b_org_y, b_y);

  Kokkos::deep_copy(h_b_x, b_x);

  typename ViewTypeA::const_type c_x = x;

  typename AT::mag_type eps = AT::epsilon() * 10;

  // Test and verify non-const input
  KokkosBlas::abs(y, x);
  Kokkos::deep_copy(h_b_y, b_y);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(h_y(i, j), AT::abs(h_x(i, j)), eps * AT::abs(h_x(i, j)));
    }
  }
  // Test and verify const input
  // Reset y
  Kokkos::deep_copy(b_y, b_org_y);
  KokkosBlas::abs(y, c_x);
  Kokkos::deep_copy(h_b_y, b_y);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(h_y(i, j), AT::abs(h_x(i, j)), eps * AT::abs(h_x(i, j)));
    }
  }
}
}  // namespace Test

template <class ScalarA, class ScalarB, class Device>
int test_abs() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_abs<view_type_a_ll, view_type_b_ll, Device>(0);
  Test::impl_test_abs<view_type_a_ll, view_type_b_ll, Device>(13);
  Test::impl_test_abs<view_type_a_ll, view_type_b_ll, Device>(1024);
  // Test::impl_test_abs<view_type_a_ll, view_type_b_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_abs<view_type_a_lr, view_type_b_lr, Device>(0);
  Test::impl_test_abs<view_type_a_lr, view_type_b_lr, Device>(13);
  Test::impl_test_abs<view_type_a_lr, view_type_b_lr, Device>(1024);
  // Test::impl_test_abs<view_type_a_lr, view_type_b_lr, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_abs<view_type_a_ls, view_type_b_ls, Device>(0);
  Test::impl_test_abs<view_type_a_ls, view_type_b_ls, Device>(13);
  Test::impl_test_abs<view_type_a_ls, view_type_b_ls, Device>(1024);
  // Test::impl_test_abs<view_type_a_ls, view_type_b_ls, Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_abs<view_type_a_ls, view_type_b_ll, Device>(1024);
  Test::impl_test_abs<view_type_a_ll, view_type_b_ls, Device>(1024);
#endif

  return 1;
}

template <class ScalarA, class ScalarB, class Device>
int test_abs_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_abs_mv<view_type_a_ll, view_type_b_ll, Device>(0, 5);
  Test::impl_test_abs_mv<view_type_a_ll, view_type_b_ll, Device>(13, 5);
  Test::impl_test_abs_mv<view_type_a_ll, view_type_b_ll, Device>(1024, 5);
  // Test::impl_test_abs_mv<view_type_a_ll, view_type_b_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_abs_mv<view_type_a_lr, view_type_b_lr, Device>(0, 5);
  Test::impl_test_abs_mv<view_type_a_lr, view_type_b_lr, Device>(13, 5);
  Test::impl_test_abs_mv<view_type_a_lr, view_type_b_lr, Device>(1024, 5);
  // Test::impl_test_abs_mv<view_type_a_lr, view_type_b_lr, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_abs_mv<view_type_a_ls, view_type_b_ls, Device>(0, 5);
  Test::impl_test_abs_mv<view_type_a_ls, view_type_b_ls, Device>(13, 5);
  Test::impl_test_abs_mv<view_type_a_ls, view_type_b_ls, Device>(1024, 5);
  // Test::impl_test_abs_mv<view_type_a_ls, view_type_b_ls, Device>(132231,5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_abs_mv<view_type_a_ls, view_type_b_ll, Device>(1024, 5);
  Test::impl_test_abs_mv<view_type_a_ll, view_type_b_ls, Device>(1024, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, abs_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_float");
  test_abs<float, float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, abs_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_mv_float");
  test_abs_mv<float, float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, abs_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_double");
  test_abs<double, double, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, abs_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_mv_double");
  test_abs_mv<double, double, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&          \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, abs_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_double");
  test_abs<Kokkos::complex<double>, Kokkos::complex<double>, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, abs_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_mv_double");
  test_abs_mv<Kokkos::complex<double>, Kokkos::complex<double>,
              TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) ||   \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, abs_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_int");
  test_abs<int, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, abs_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_mv_int");
  test_abs_mv<int, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

/*#if !defined(KOKKOSKERNELS_ETI_ONLY) &&
!defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS) TEST_F( TestCategory,
abs_double_int ) { test_abs<double,int,TestExecSpace> ();
}
TEST_F( TestCategory, abs_double_mv_int ) {
    test_abs_mv<double,int,TestExecSpace> ();
}
#endif*/
