#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas1_nrm2w.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class Device>
void impl_test_nrm2w(int N) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef Kokkos::ArithTraits<ScalarA> AT;

  ViewTypeA a("A", N);
  ViewTypeA w("W", N);

  typename ViewTypeA::HostMirror h_a = Kokkos::create_mirror_view(a);
  typename ViewTypeA::HostMirror h_w = Kokkos::create_mirror_view(w);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  ScalarA randStart, randEnd;
  Test::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(a, rand_pool, randStart, randEnd);
  Kokkos::fill_random(w, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(h_a, a);
  Kokkos::deep_copy(h_w, w);

  double eps = std::is_same<ScalarA, float>::value ? 2 * 1e-5 : 1e-7;

  typename AT::mag_type expected_result = 0;
  for (int i = 0; i < N; i++) {
    typename AT::mag_type term = AT::abs(h_a(i)) / AT::abs(h_w(i));
    expected_result += term * term;
  }
  expected_result =
      Kokkos::ArithTraits<typename AT::mag_type>::sqrt(expected_result);

  typename AT::mag_type nonconst_result = KokkosBlas::nrm2w(a, w);
  EXPECT_NEAR_KK(nonconst_result, expected_result, eps * expected_result);
}

template <class ViewTypeA, class Device>
void impl_test_nrm2w_mv(int N, int K) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef Kokkos::ArithTraits<ScalarA> AT;

  typedef multivector_layout_adapter<ViewTypeA> vfA_type;

  typename vfA_type::BaseType b_a("A", N, K);
  typename vfA_type::BaseType b_w("W", N, K);

  ViewTypeA a = vfA_type::view(b_a);
  ViewTypeA w = vfA_type::view(b_w);

  typedef multivector_layout_adapter<typename ViewTypeA::HostMirror> h_vfA_type;

  typename h_vfA_type::BaseType h_b_a = Kokkos::create_mirror_view(b_a);
  typename h_vfA_type::BaseType h_b_w = Kokkos::create_mirror_view(b_w);

  typename ViewTypeA::HostMirror h_a = h_vfA_type::view(h_b_a);
  typename ViewTypeA::HostMirror h_w = h_vfA_type::view(h_b_w);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  ScalarA randStart, randEnd;
  Test::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(b_a, rand_pool, randStart, randEnd);
  Kokkos::fill_random(b_w, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(h_b_a, b_a);
  Kokkos::deep_copy(h_b_w, b_w);

  typename AT::mag_type* expected_result = new typename AT::mag_type[K];
  for (int j = 0; j < K; j++) {
    expected_result[j] = typename AT::mag_type();
    for (int i = 0; i < N; i++) {
      typename AT::mag_type term = AT::abs(h_a(i, j)) / AT::abs(h_w(i, j));
      expected_result[j] += term * term;
    }
    expected_result[j] =
        Kokkos::ArithTraits<typename AT::mag_type>::sqrt(expected_result[j]);
  }

  double eps = std::is_same<ScalarA, float>::value ? 2 * 1e-5 : 1e-7;

  Kokkos::View<typename AT::mag_type*, Device> r("Dot::Result", K);
  KokkosBlas::nrm2w(r, a, w);
  auto r_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), r);

  for (int k = 0; k < K; k++) {
    typename AT::mag_type nonconst_result = r_host(k);
    EXPECT_NEAR_KK(nonconst_result, expected_result[k],
                   eps * expected_result[k]);
  }

  delete[] expected_result;
}
}  // namespace Test

template <class ScalarA, class Device>
int test_nrm2w() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrm2w<view_type_a_ll, Device>(0);
  Test::impl_test_nrm2w<view_type_a_ll, Device>(13);
  Test::impl_test_nrm2w<view_type_a_ll, Device>(1024);
  // Test::impl_test_nrm2<view_type_a_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrm2w<view_type_a_lr, Device>(0);
  Test::impl_test_nrm2w<view_type_a_lr, Device>(13);
  Test::impl_test_nrm2w<view_type_a_lr, Device>(1024);
  // Test::impl_test_nrm2<view_type_a_lr, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_nrm2w<view_type_a_ls, Device>(0);
  Test::impl_test_nrm2w<view_type_a_ls, Device>(13);
  Test::impl_test_nrm2w<view_type_a_ls, Device>(1024);
  // Test::impl_test_nrm2<view_type_a_ls, Device>(132231);
#endif

  return 1;
}

template <class ScalarA, class Device>
int test_nrm2w_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrm2w_mv<view_type_a_ll, Device>(0, 5);
  Test::impl_test_nrm2w_mv<view_type_a_ll, Device>(13, 5);
  Test::impl_test_nrm2w_mv<view_type_a_ll, Device>(1024, 5);
  Test::impl_test_nrm2w_mv<view_type_a_ll, Device>(789, 1);
  // Test::impl_test_nrm2w_mv<view_type_a_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrm2w_mv<view_type_a_lr, Device>(0, 5);
  Test::impl_test_nrm2w_mv<view_type_a_lr, Device>(13, 5);
  Test::impl_test_nrm2w_mv<view_type_a_lr, Device>(1024, 5);
  Test::impl_test_nrm2w_mv<view_type_a_lr, Device>(789, 1);
  // Test::impl_test_nrm2w_mv<view_type_a_lr, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_nrm2w_mv<view_type_a_ls, Device>(0, 5);
  Test::impl_test_nrm2w_mv<view_type_a_ls, Device>(13, 5);
  Test::impl_test_nrm2w_mv<view_type_a_ls, Device>(1024, 5);
  Test::impl_test_nrm2w_mv<view_type_a_ls, Device>(789, 1);
  // Test::impl_test_nrm2w_mv<view_type_a_ls, Device>(132231,5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm2w_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_float");
  test_nrm2w<float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm2w_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_mv_float");
  test_nrm2w_mv<float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm2w_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_double");
  test_nrm2w<double, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm2w_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_mv_double");
  test_nrm2w_mv<double, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&          \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm2w_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_complex_double");
  test_nrm2w<Kokkos::complex<double>, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm2w_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_mv_complex_double");
  test_nrm2w_mv<Kokkos::complex<double>, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) ||   \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm2w_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_int");
  test_nrm2w<int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm2w_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_mv_int");
  test_nrm2w_mv<int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif
