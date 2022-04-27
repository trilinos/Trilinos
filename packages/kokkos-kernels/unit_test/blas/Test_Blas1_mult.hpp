#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas1_mult.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class Device>
void impl_test_mult(int N) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef typename ViewTypeC::value_type ScalarC;

  ScalarA a  = 3;
  ScalarB b  = 5;
  double eps = std::is_same<ScalarC, float>::value ? 1e-4 : 1e-7;

  ViewTypeA x("X", N);
  ViewTypeB y("Y", N);
  ViewTypeC z("Y", N);
  ViewTypeC b_org_z("Org_Z", N);

  typename ViewTypeA::const_type c_x = x;
  typename ViewTypeB::const_type c_y = y;

  typename ViewTypeA::HostMirror h_x = Kokkos::create_mirror_view(x);
  typename ViewTypeB::HostMirror h_y = Kokkos::create_mirror_view(y);
  typename ViewTypeC::HostMirror h_z = Kokkos::create_mirror_view(z);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(x, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(y, rand_pool, randStart, randEnd);
  }
  {
    ScalarC randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(z, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(b_org_z, z);
  auto h_b_org_z =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b_org_z);

  Kokkos::deep_copy(h_x, x);
  Kokkos::deep_copy(h_y, y);

  // expected_result = ScalarC(b*h_z(i) + a*h_x(i)*h_y(i))

  KokkosBlas::mult(b, z, a, x, y);
  Kokkos::deep_copy(h_z, z);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(a * h_x(i) * h_y(i) + b * h_b_org_z(i), h_z(i), eps);
  }

  Kokkos::deep_copy(z, b_org_z);
  KokkosBlas::mult(b, z, a, x, c_y);
  Kokkos::deep_copy(h_z, z);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(a * h_x(i) * h_y(i) + b * h_b_org_z(i), h_z(i), eps);
  }

  Kokkos::deep_copy(z, b_org_z);
  KokkosBlas::mult(b, z, a, c_x, c_y);
  Kokkos::deep_copy(h_z, z);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(a * h_x(i) * h_y(i) + b * h_b_org_z(i), h_z(i), eps);
  }
}

template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class Device>
void impl_test_mult_mv(int N, int K) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef typename ViewTypeC::value_type ScalarC;

  typedef multivector_layout_adapter<ViewTypeB> vfB_type;
  typedef multivector_layout_adapter<ViewTypeC> vfC_type;

  ViewTypeA x("X", N);
  typename vfB_type::BaseType b_y("Y", N, K);
  typename vfC_type::BaseType b_z("Z", N, K);
  typename vfC_type::BaseType b_org_z("Z", N, K);

  ViewTypeB y = vfB_type::view(b_y);
  ViewTypeC z = vfC_type::view(b_z);

  typedef multivector_layout_adapter<typename ViewTypeB::HostMirror> h_vfB_type;
  typedef multivector_layout_adapter<typename ViewTypeC::HostMirror> h_vfC_type;

  typename h_vfB_type::BaseType h_b_y = Kokkos::create_mirror_view(b_y);
  typename h_vfC_type::BaseType h_b_z = Kokkos::create_mirror_view(b_z);

  typename ViewTypeA::HostMirror h_x = Kokkos::create_mirror_view(x);
  typename ViewTypeB::HostMirror h_y = h_vfB_type::view(h_b_y);
  typename ViewTypeC::HostMirror h_z = h_vfC_type::view(h_b_z);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(x, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(b_y, rand_pool, randStart, randEnd);
  }
  {
    ScalarC randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(b_z, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(b_org_z, b_z);
  auto h_b_org_z =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b_org_z);

  Kokkos::deep_copy(h_x, x);
  Kokkos::deep_copy(h_b_y, b_y);
  Kokkos::deep_copy(h_b_z, b_z);

  ScalarA a                          = 3;
  ScalarB b                          = 5;
  typename ViewTypeA::const_type c_x = x;
  typename ViewTypeB::const_type c_y = y;

  double eps = std::is_same<ScalarA, float>::value ? 1e-4 : 1e-7;

  KokkosBlas::mult(b, z, a, x, y);
  Kokkos::deep_copy(h_b_z, b_z);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(a * h_x(i) * h_y(i, j) + b * h_b_org_z(i, j), h_z(i, j),
                     eps);
    }
  }

  Kokkos::deep_copy(b_z, b_org_z);
  KokkosBlas::mult(b, z, a, x, c_y);
  Kokkos::deep_copy(h_b_z, b_z);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(a * h_x(i) * h_y(i, j) + b * h_b_org_z(i, j), h_z(i, j),
                     eps);
    }
  }
}
}  // namespace Test

template <class ScalarA, class ScalarB, class ScalarC, class Device>
int test_mult() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  typedef Kokkos::View<ScalarC*, Kokkos::LayoutLeft, Device> view_type_c_ll;
  Test::impl_test_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(
      0);
  Test::impl_test_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(
      13);
  Test::impl_test_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(
      1024);
  // Test::impl_test_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll,
  // Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutRight, Device> view_type_b_lr;
  typedef Kokkos::View<ScalarC*, Kokkos::LayoutRight, Device> view_type_c_lr;
  Test::impl_test_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(
      0);
  Test::impl_test_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(
      13);
  Test::impl_test_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(
      1024);
  // Test::impl_test_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr,
  // Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, Device> view_type_b_ls;
  typedef Kokkos::View<ScalarC*, Kokkos::LayoutStride, Device> view_type_c_ls;
  Test::impl_test_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(
      0);
  Test::impl_test_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(
      13);
  Test::impl_test_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(
      1024);
  // Test::impl_test_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls,
  // Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_mult<view_type_a_ls, view_type_b_ll, view_type_c_lr, Device>(
      1024);
  Test::impl_test_mult<view_type_a_ll, view_type_b_ls, view_type_c_lr, Device>(
      1024);
#endif

  return 1;
}

template <class ScalarA, class ScalarB, class ScalarC, class Device>
int test_mult_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  typedef Kokkos::View<ScalarC**, Kokkos::LayoutLeft, Device> view_type_c_ll;
  Test::impl_test_mult_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll,
                          Device>(0, 5);
  Test::impl_test_mult_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll,
                          Device>(13, 5);
  Test::impl_test_mult_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll,
                          Device>(1024, 5);
  // Test::impl_test_mult_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll,
  // Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  typedef Kokkos::View<ScalarC**, Kokkos::LayoutRight, Device> view_type_c_lr;
  Test::impl_test_mult_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr,
                          Device>(0, 5);
  Test::impl_test_mult_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr,
                          Device>(13, 5);
  Test::impl_test_mult_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr,
                          Device>(1024, 5);
  // Test::impl_test_mult_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr,
  // Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutStride, Device> view_type_b_ls;
  typedef Kokkos::View<ScalarC**, Kokkos::LayoutStride, Device> view_type_c_ls;
  Test::impl_test_mult_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls,
                          Device>(0, 5);
  Test::impl_test_mult_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls,
                          Device>(13, 5);
  Test::impl_test_mult_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls,
                          Device>(1024, 5);
  // Test::impl_test_mult_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls,
  // Device>(132231,5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_mult_mv<view_type_a_ls, view_type_b_ll, view_type_c_lr,
                          Device>(1024, 5);
  Test::impl_test_mult_mv<view_type_a_ll, view_type_b_ls, view_type_c_lr,
                          Device>(1024, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, mult_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::mult_float");
  test_mult<float, float, float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, mult_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::mult_float");
  test_mult_mv<float, float, float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, mult_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::mult_double");
  test_mult<double, double, double, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, mult_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::mult_mv_double");
  test_mult_mv<double, double, double, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&          \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, mult_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::mult_complex_double");
  test_mult<Kokkos::complex<double>, Kokkos::complex<double>,
            Kokkos::complex<double>, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, mult_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::mult_mv_complex_double");
  test_mult_mv<Kokkos::complex<double>, Kokkos::complex<double>,
               Kokkos::complex<double>, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) ||   \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, mult_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::mult_int");
  test_mult<int, int, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, mult_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::mult_mv_int");
  test_mult_mv<int, int, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, mult_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::mult_double_int");
  test_mult<double, int, float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, mult_mv_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::mult_mv_double_int");
  test_mult_mv<double, int, float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif
