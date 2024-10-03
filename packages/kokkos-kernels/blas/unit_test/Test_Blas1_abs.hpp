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
#include <KokkosBlas1_abs.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_abs(int N) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef Kokkos::ArithTraits<ScalarA> AT;

  typename AT::mag_type eps = AT::epsilon() * 10;

  view_stride_adapter<ViewTypeA> x("X", N);
  view_stride_adapter<ViewTypeB> y("Y", N);
  view_stride_adapter<ViewTypeB> org_y("Org_Y", N);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(y.d_view, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(org_y.h_base, y.d_base);

  Kokkos::deep_copy(x.h_base, x.d_base);

  // Run with nonconst input
  KokkosBlas::abs(y.d_view, x.d_view);
  // Copy result to host (h_y is subview of h_b_y)
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(y.h_view(i), AT::abs(x.h_view(i)), eps * AT::abs(x.h_view(i)));
  }
  // Run with const input
  // Reset output
  Kokkos::deep_copy(y.d_base, org_y.h_base);
  KokkosBlas::abs(y.d_view, x.d_view_const);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(y.h_view(i), AT::abs(x.h_view(i)), eps * AT::abs(x.h_view(i)));
  }
}

template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_abs_mv(int N, int K) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef Kokkos::ArithTraits<ScalarA> AT;

  view_stride_adapter<ViewTypeA> x("X", N, K);
  view_stride_adapter<ViewTypeB> y("Y", N, K);
  view_stride_adapter<ViewTypeB> org_y("Org_Y", N, K);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(y.d_view, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(org_y.h_base, y.d_base);

  Kokkos::deep_copy(x.h_base, x.d_base);

  typename AT::mag_type eps = AT::epsilon() * 10;

  // Test and verify non-const input
  KokkosBlas::abs(y.d_view, x.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(y.h_view(i, j), AT::abs(x.h_view(i, j)), eps * AT::abs(x.h_view(i, j)));
    }
  }
  // Test and verify const input
  // Reset y
  Kokkos::deep_copy(y.d_base, org_y.h_base);
  KokkosBlas::abs(y.d_view, x.d_view_const);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(y.h_view(i, j), AT::abs(x.h_view(i, j)), eps * AT::abs(x.h_view(i, j)));
    }
  }
}
}  // namespace Test

template <class ScalarA, class ScalarB, class Device>
int test_abs() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_abs<view_type_a_ll, view_type_b_ll, Device>(0);
  Test::impl_test_abs<view_type_a_ll, view_type_b_ll, Device>(13);
  Test::impl_test_abs<view_type_a_ll, view_type_b_ll, Device>(1024);
  // Test::impl_test_abs<view_type_a_ll, view_type_b_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_abs<view_type_a_lr, view_type_b_lr, Device>(0);
  Test::impl_test_abs<view_type_a_lr, view_type_b_lr, Device>(13);
  Test::impl_test_abs<view_type_a_lr, view_type_b_lr, Device>(1024);
  // Test::impl_test_abs<view_type_a_lr, view_type_b_lr, Device>(132231);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_abs<view_type_a_ls, view_type_b_ls, Device>(0);
  Test::impl_test_abs<view_type_a_ls, view_type_b_ls, Device>(13);
  Test::impl_test_abs<view_type_a_ls, view_type_b_ls, Device>(1024);
  // Test::impl_test_abs<view_type_a_ls, view_type_b_ls, Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_abs<view_type_a_ls, view_type_b_ll, Device>(1024);
  Test::impl_test_abs<view_type_a_ll, view_type_b_ls, Device>(1024);
#endif

  return 1;
}

template <class ScalarA, class ScalarB, class Device>
int test_abs_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_abs_mv<view_type_a_ll, view_type_b_ll, Device>(0, 5);
  Test::impl_test_abs_mv<view_type_a_ll, view_type_b_ll, Device>(13, 5);
  Test::impl_test_abs_mv<view_type_a_ll, view_type_b_ll, Device>(1024, 5);
  // Test::impl_test_abs_mv<view_type_a_ll, view_type_b_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_abs_mv<view_type_a_lr, view_type_b_lr, Device>(0, 5);
  Test::impl_test_abs_mv<view_type_a_lr, view_type_b_lr, Device>(13, 5);
  Test::impl_test_abs_mv<view_type_a_lr, view_type_b_lr, Device>(1024, 5);
  // Test::impl_test_abs_mv<view_type_a_lr, view_type_b_lr, Device>(132231,5);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_abs_mv<view_type_a_ls, view_type_b_ls, Device>(0, 5);
  Test::impl_test_abs_mv<view_type_a_ls, view_type_b_ls, Device>(13, 5);
  Test::impl_test_abs_mv<view_type_a_ls, view_type_b_ls, Device>(1024, 5);
  // Test::impl_test_abs_mv<view_type_a_ls, view_type_b_ls, Device>(132231,5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_abs_mv<view_type_a_ls, view_type_b_ll, Device>(1024, 5);
  Test::impl_test_abs_mv<view_type_a_ll, view_type_b_ls, Device>(1024, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, abs_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_float");
  test_abs<float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, abs_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_mv_float");
  test_abs_mv<float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, abs_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_double");
  test_abs<double, double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, abs_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_mv_double");
  test_abs_mv<double, double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, abs_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_double");
  test_abs<Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, abs_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_mv_double");
  test_abs_mv<Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, abs_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_int");
  test_abs<int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, abs_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::abs_mv_int");
  test_abs_mv<int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

/*#if !defined(KOKKOSKERNELS_ETI_ONLY) &&
!defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS) TEST_F( TestCategory,
abs_double_int ) { test_abs<double,int,TestDevice> ();
}
TEST_F( TestCategory, abs_double_mv_int ) {
    test_abs_mv<double,int,TestDevice> ();
}
#endif*/
