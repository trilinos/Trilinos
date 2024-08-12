//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class Device>
void impl_test_nrm1(int N) {
  using ScalarA  = typename ViewTypeA::value_type;
  using AT       = Kokkos::ArithTraits<ScalarA>;
  using mag_type = typename AT::mag_type;
  using MAT      = Kokkos::ArithTraits<mag_type>;

  view_stride_adapter<ViewTypeA> a("a", N);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  ScalarA randStart, randEnd;
  Test::getRandomBounds(10.0, randStart, randEnd);
  Kokkos::fill_random(a.d_view, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(a.h_base, a.d_base);

  double eps = (std::is_same<typename Kokkos::ArithTraits<ScalarA>::mag_type, float>::value ? 1e-4 : 1e-7);

  mag_type expected_result = 0;
  for (int i = 0; i < N; i++) {
    // note: for complex, BLAS asum (aka our nrm1) is _not_
    // the sum of magnitudes - it's the sum of absolute real and imaginary
    // parts. See netlib, MKL, and CUBLAS documentation.
    //
    // This is safe; ArithTraits<T>::imag is 0 if T is real.
    expected_result += MAT::abs(AT::real(a.h_view(i))) + MAT::abs(AT::imag(a.h_view(i)));
  }

  mag_type nonconst_result = KokkosBlas::nrm1(a.d_view);
  EXPECT_NEAR_KK(nonconst_result, expected_result, eps * expected_result);

  mag_type const_result = KokkosBlas::nrm1(a.d_view_const);
  EXPECT_NEAR_KK(const_result, expected_result, eps * expected_result);
}

template <class ViewTypeA, class Device>
void impl_test_nrm1_mv(int N, int K) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef Kokkos::ArithTraits<ScalarA> AT;
  typedef typename AT::mag_type mag_type;
  typedef Kokkos::ArithTraits<mag_type> MAT;

  view_stride_adapter<ViewTypeA> a("A", N, K);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  ScalarA randStart, randEnd;
  Test::getRandomBounds(10.0, randStart, randEnd);
  Kokkos::fill_random(a.d_view, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(a.h_base, a.d_base);

  double eps = (std::is_same<typename Kokkos::ArithTraits<ScalarA>::mag_type, float>::value ? 1e-4 : 1e-7);

  Kokkos::View<mag_type*, Kokkos::HostSpace> expected_result("Expected Nrm1", K);
  for (int k = 0; k < K; k++) {
    expected_result(k) = MAT::zero();
    for (int i = 0; i < N; i++) {
      expected_result(k) += MAT::abs(AT::real(a.h_view(i, k))) + MAT::abs(AT::imag(a.h_view(i, k)));
    }
  }

  Kokkos::View<mag_type*, Kokkos::HostSpace> r("Nrm1::Result", K);
  Kokkos::View<mag_type*, Kokkos::HostSpace> c_r("Nrm1::ConstResult", K);

  KokkosBlas::nrm1(r, a.d_view);
  KokkosBlas::nrm1(c_r, a.d_view_const);
  Kokkos::fence();
  for (int k = 0; k < K; k++) {
    EXPECT_NEAR_KK(r(k), expected_result(k), eps * expected_result(k));
  }
}
}  // namespace Test

template <class ScalarA, class Device>
int test_nrm1() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrm1<view_type_a_ll, Device>(0);
  Test::impl_test_nrm1<view_type_a_ll, Device>(13);
  Test::impl_test_nrm1<view_type_a_ll, Device>(1024);
  Test::impl_test_nrm1<view_type_a_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrm1<view_type_a_lr, Device>(0);
  Test::impl_test_nrm1<view_type_a_lr, Device>(13);
  Test::impl_test_nrm1<view_type_a_lr, Device>(1024);
  Test::impl_test_nrm1<view_type_a_lr, Device>(132231);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
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
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrm1_mv<view_type_a_ll, Device>(0, 5);
  Test::impl_test_nrm1_mv<view_type_a_ll, Device>(13, 5);
  Test::impl_test_nrm1_mv<view_type_a_ll, Device>(1024, 5);
  Test::impl_test_nrm1_mv<view_type_a_ll, Device>(789, 1);
  Test::impl_test_nrm1_mv<view_type_a_ll, Device>(132231, 5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrm1_mv<view_type_a_lr, Device>(0, 5);
  Test::impl_test_nrm1_mv<view_type_a_lr, Device>(13, 5);
  Test::impl_test_nrm1_mv<view_type_a_lr, Device>(1024, 5);
  Test::impl_test_nrm1_mv<view_type_a_lr, Device>(789, 1);
  Test::impl_test_nrm1_mv<view_type_a_lr, Device>(132231, 5);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
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
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm1_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_float");
  test_nrm1<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm1_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_mv_float");
  test_nrm1_mv<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm1_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_double");
  test_nrm1<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm1_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_mv_double");
  test_nrm1_mv<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm1_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_complex_double");
  test_nrm1<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm1_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_mv_complex_double");
  test_nrm1_mv<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm1_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_int");
  test_nrm1<int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm1_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm1_mv_int");
  test_nrm1_mv<int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
