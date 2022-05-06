#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class Device>
void impl_test_nrm1(int N) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef Kokkos::ArithTraits<ScalarA> AT;
  typedef typename AT::mag_type mag_type;
  typedef Kokkos::ArithTraits<mag_type> MAT;

  ViewTypeA a("A", N);

  typename ViewTypeA::HostMirror h_a = Kokkos::create_mirror_view(a);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  ScalarA randStart, randEnd;
  Test::getRandomBounds(10.0, randStart, randEnd);
  Kokkos::fill_random(a, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(h_a, a);

  typename ViewTypeA::const_type c_a = a;
  double eps = (std::is_same<typename Kokkos::ArithTraits<ScalarA>::mag_type,
                             float>::value
                    ? 1e-4
                    : 1e-7);

  mag_type expected_result = 0;
  for (int i = 0; i < N; i++) {
    // note: for complex, BLAS asum (aka our nrm1) is _not_
    // the sum of magnitudes - it's the sum of absolute real and imaginary
    // parts. See netlib, MKL, and CUBLAS documentation.
    //
    // This is safe; ArithTraits<T>::imag is 0 if T is real.
    expected_result += MAT::abs(AT::real(h_a(i))) + MAT::abs(AT::imag(h_a(i)));
  }

  mag_type nonconst_result = KokkosBlas::nrm1(a);
  EXPECT_NEAR_KK(nonconst_result, expected_result, eps * expected_result);

  mag_type const_result = KokkosBlas::nrm1(c_a);
  EXPECT_NEAR_KK(const_result, expected_result, eps * expected_result);
}

template <class ViewTypeA, class Device>
void impl_test_nrm1_mv(int N, int K) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef Kokkos::Details::ArithTraits<ScalarA> AT;
  typedef typename AT::mag_type mag_type;
  typedef Kokkos::ArithTraits<mag_type> MAT;

  typedef multivector_layout_adapter<ViewTypeA> vfA_type;

  typename vfA_type::BaseType b_a("A", N, K);

  ViewTypeA a = vfA_type::view(b_a);

  typedef multivector_layout_adapter<typename ViewTypeA::HostMirror> h_vfA_type;

  typename h_vfA_type::BaseType h_b_a = Kokkos::create_mirror_view(b_a);

  typename ViewTypeA::HostMirror h_a = h_vfA_type::view(h_b_a);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  ScalarA randStart, randEnd;
  Test::getRandomBounds(10.0, randStart, randEnd);
  Kokkos::fill_random(b_a, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(h_b_a, b_a);

  typename ViewTypeA::const_type c_a = a;

  double eps = (std::is_same<typename Kokkos::ArithTraits<ScalarA>::mag_type,
                             float>::value
                    ? 1e-4
                    : 1e-7);

  Kokkos::View<mag_type*, Kokkos::HostSpace> expected_result("Expected Nrm1",
                                                             K);
  for (int k = 0; k < K; k++) {
    expected_result(k) = MAT::zero();
    for (int i = 0; i < N; i++) {
      expected_result(k) +=
          MAT::abs(AT::real(h_a(i, k))) + MAT::abs(AT::imag(h_a(i, k)));
    }
  }

  Kokkos::View<mag_type*, Kokkos::HostSpace> r("Nrm1::Result", K);
  Kokkos::View<mag_type*, Kokkos::HostSpace> c_r("Nrm1::ConstResult", K);

  KokkosBlas::nrm1(r, a);
  KokkosBlas::nrm1(c_r, a);
  for (int k = 0; k < K; k++) {
    EXPECT_NEAR_KK(r(k), expected_result(k), eps * expected_result(k));
    EXPECT_NEAR_KK(c_r(k), expected_result(k), eps * expected_result(k));
  }
}
}  // namespace Test

template <class ScalarA, class Device>
int test_nrm1() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrm1<view_type_a_ll, Device>(0);
  Test::impl_test_nrm1<view_type_a_ll, Device>(13);
  Test::impl_test_nrm1<view_type_a_ll, Device>(1024);
  Test::impl_test_nrm1<view_type_a_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrm1<view_type_a_lr, Device>(0);
  Test::impl_test_nrm1<view_type_a_lr, Device>(13);
  Test::impl_test_nrm1<view_type_a_lr, Device>(1024);
  Test::impl_test_nrm1<view_type_a_lr, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_nrm1<view_type_a_ls, Device>(0);
  Test::impl_test_nrm1<view_type_a_ls, Device>(13);
  Test::impl_test_nrm1<view_type_a_ls, Device>(1024);
  Test::impl_test_nrm1<view_type_a_ls, Device>(132231);
#endif

  return 1;
}

template <class ScalarA, class Device>
int test_nrm1_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrm1_mv<view_type_a_ll, Device>(0, 5);
  Test::impl_test_nrm1_mv<view_type_a_ll, Device>(13, 5);
  Test::impl_test_nrm1_mv<view_type_a_ll, Device>(1024, 5);
  Test::impl_test_nrm1_mv<view_type_a_ll, Device>(789, 1);
  Test::impl_test_nrm1_mv<view_type_a_ll, Device>(132231, 5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrm1_mv<view_type_a_lr, Device>(0, 5);
  Test::impl_test_nrm1_mv<view_type_a_lr, Device>(13, 5);
  Test::impl_test_nrm1_mv<view_type_a_lr, Device>(1024, 5);
  Test::impl_test_nrm1_mv<view_type_a_lr, Device>(789, 1);
  Test::impl_test_nrm1_mv<view_type_a_lr, Device>(132231, 5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_nrm1_mv<view_type_a_ls, Device>(0, 5);
  Test::impl_test_nrm1_mv<view_type_a_ls, Device>(13, 5);
  Test::impl_test_nrm1_mv<view_type_a_ls, Device>(1024, 5);
  Test::impl_test_nrm1_mv<view_type_a_ls, Device>(789, 1);
  Test::impl_test_nrm1_mv<view_type_a_ls, Device>(132231, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm1_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_float");
  test_nrm1<float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm1_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_mv_float");
  test_nrm1_mv<float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm1_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_double");
  test_nrm1<double, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm1_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_mv_double");
  test_nrm1_mv<double, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&          \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm1_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_complex_double");
  test_nrm1<Kokkos::complex<double>, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm1_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_mv_complex_double");
  test_nrm1_mv<Kokkos::complex<double>, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) ||   \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm1_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_int");
  test_nrm1<int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm1_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_mv_int");
  test_nrm1_mv<int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif
