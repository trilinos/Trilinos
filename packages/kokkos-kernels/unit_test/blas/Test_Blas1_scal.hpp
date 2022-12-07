#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas1_scal.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_scal(int N) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef Kokkos::Details::ArithTraits<ScalarA> AT;

  ScalarA a(3);
  typename AT::mag_type eps = AT::epsilon() * 1000;

  ViewTypeA x("X", N);
  ViewTypeB y("Y", N);
  ViewTypeB org_y("Org_Y", N);

  typename ViewTypeA::const_type c_x = x;
  typename ViewTypeB::const_type c_y = y;

  typename ViewTypeA::HostMirror h_x = Kokkos::create_mirror_view(x);
  typename ViewTypeB::HostMirror h_y = Kokkos::create_mirror_view(y);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(x, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(y, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(org_y, y);

  Kokkos::deep_copy(h_x, x);

  KokkosBlas::scal(y, a, x);
  Kokkos::deep_copy(h_y, y);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(a * h_x(i), h_y(i), eps);
  }

  Kokkos::deep_copy(y, org_y);
  KokkosBlas::scal(y, a, c_x);
  Kokkos::deep_copy(h_y, y);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(a * h_x(i), h_y(i), eps);
  }
}

template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_scal_mv(int N, int K) {
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

  Kokkos::fence();

  Kokkos::deep_copy(b_org_y, b_y);

  Kokkos::deep_copy(h_b_x, b_x);

  ScalarA a(3.0);
  typename ViewTypeA::const_type c_x = x;

  typename AT::mag_type eps = AT::epsilon() * 1000;

  Kokkos::View<ScalarB*, Kokkos::HostSpace> r("Dot::Result", K);

  KokkosBlas::scal(y, a, x);
  Kokkos::deep_copy(h_b_y, b_y);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(a * h_x(i, j), h_y(i, j), eps);
    }
  }

  Kokkos::deep_copy(b_y, b_org_y);
  KokkosBlas::scal(y, a, c_x);
  Kokkos::deep_copy(h_b_y, b_y);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(a * h_x(i, j), h_y(i, j), eps);
    }
  }

  // Generate 'params' view with dimension == number of multivectors; each entry
  // will be different scalar to scale y
  Kokkos::View<ScalarA*, Device> params("Params", K);
  for (int j = 0; j < K; j++) {
    Kokkos::View<ScalarA, Device> param_j(params, j);
    Kokkos::deep_copy(param_j, ScalarA(3 + j));
  }

  auto h_params =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), params);

  KokkosBlas::scal(y, params, x);
  Kokkos::deep_copy(h_b_y, b_y);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(h_params(j) * h_x(i, j), h_y(i, j), eps);
    }
  }

  Kokkos::deep_copy(b_y, b_org_y);
  KokkosBlas::scal(y, params, c_x);
  Kokkos::deep_copy(h_b_y, b_y);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(h_params(j) * h_x(i, j), h_y(i, j), eps);
    }
  }
}
}  // namespace Test

template <class ScalarA, class ScalarB, class Device>
int test_scal() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_scal<view_type_a_ll, view_type_b_ll, Device>(0);
  Test::impl_test_scal<view_type_a_ll, view_type_b_ll, Device>(13);
  Test::impl_test_scal<view_type_a_ll, view_type_b_ll, Device>(1024);
  // Test::impl_test_scal<view_type_a_ll, view_type_b_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_scal<view_type_a_lr, view_type_b_lr, Device>(0);
  Test::impl_test_scal<view_type_a_lr, view_type_b_lr, Device>(13);
  Test::impl_test_scal<view_type_a_lr, view_type_b_lr, Device>(1024);
  // Test::impl_test_scal<view_type_a_lr, view_type_b_lr, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_scal<view_type_a_ls, view_type_b_ls, Device>(0);
  Test::impl_test_scal<view_type_a_ls, view_type_b_ls, Device>(13);
  Test::impl_test_scal<view_type_a_ls, view_type_b_ls, Device>(1024);
  // Test::impl_test_scal<view_type_a_ls, view_type_b_ls, Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_scal<view_type_a_ls, view_type_b_ll, Device>(1024);
  Test::impl_test_scal<view_type_a_ll, view_type_b_ls, Device>(1024);
#endif

  return 1;
}

template <class ScalarA, class ScalarB, class Device>
int test_scal_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_scal_mv<view_type_a_ll, view_type_b_ll, Device>(0, 5);
  Test::impl_test_scal_mv<view_type_a_ll, view_type_b_ll, Device>(13, 5);
  Test::impl_test_scal_mv<view_type_a_ll, view_type_b_ll, Device>(1024, 5);
  // Test::impl_test_scal_mv<view_type_a_ll, view_type_b_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_scal_mv<view_type_a_lr, view_type_b_lr, Device>(0, 5);
  Test::impl_test_scal_mv<view_type_a_lr, view_type_b_lr, Device>(13, 5);
  Test::impl_test_scal_mv<view_type_a_lr, view_type_b_lr, Device>(1024, 5);
  // Test::impl_test_scal_mv<view_type_a_lr, view_type_b_lr, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_scal_mv<view_type_a_ls, view_type_b_ls, Device>(0, 5);
  Test::impl_test_scal_mv<view_type_a_ls, view_type_b_ls, Device>(13, 5);
  Test::impl_test_scal_mv<view_type_a_ls, view_type_b_ls, Device>(1024, 5);
  // Test::impl_test_scal_mv<view_type_a_ls, view_type_b_ls, Device>(132231,5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_scal_mv<view_type_a_ls, view_type_b_ll, Device>(1024, 5);
  Test::impl_test_scal_mv<view_type_a_ll, view_type_b_ls, Device>(1024, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, scal_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_float");
  test_scal<float, float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, scal_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_mv_float");
  test_scal_mv<float, float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, scal_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_double");
  test_scal<double, double, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, scal_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_mv_double");
  test_scal_mv<double, double, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&          \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, scal_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_complex_double");
  test_scal<Kokkos::complex<double>, Kokkos::complex<double>, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, scal_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_mv_complex_double");
  test_scal_mv<Kokkos::complex<double>, Kokkos::complex<double>,
               TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) ||   \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, scal_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_int");
  test_scal<int, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, scal_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_mv_int");
  test_scal_mv<int, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, scal_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_double_int");
  test_scal<double, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, scal_mv_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_mv_double_int");
  test_scal_mv<double, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif
