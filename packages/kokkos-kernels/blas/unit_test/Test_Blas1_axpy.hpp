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
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_axpy(int N) {
  using ScalarA    = typename ViewTypeA::value_type;
  using ScalarB    = typename ViewTypeB::value_type;
  using MagnitudeB = typename Kokkos::ArithTraits<ScalarB>::mag_type;

  ScalarA a                = 3;
  const MagnitudeB max_val = 10;
  const MagnitudeB eps     = Kokkos::ArithTraits<ScalarB>::epsilon();
  const MagnitudeB max_error =
      (static_cast<MagnitudeB>(Kokkos::ArithTraits<ScalarA>::abs(a)) * max_val + max_val) * eps;

  view_stride_adapter<ViewTypeA> x("X", N);
  view_stride_adapter<ViewTypeB> y("Y", N);
  view_stride_adapter<ViewTypeB> org_y("Org_Y", N);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(max_val, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(max_val, randStart, randEnd);
    Kokkos::fill_random(y.d_view, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(x.h_base, x.d_base);
  Kokkos::deep_copy(org_y.h_base, y.d_base);

  KokkosBlas::axpy(a, x.d_view, y.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);

  for (int i = 0; i < N; i++) {
    ScalarB expected = a * x.h_view(i) + org_y.h_view(i);
    EXPECT_NEAR_KK(expected, y.h_view(i), 2 * max_error);
  }

  // reset y to orig, and run again with const-valued x
  Kokkos::deep_copy(y.d_base, org_y.h_base);
  KokkosBlas::axpy(a, x.d_view_const, y.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    ScalarB expected = a * x.h_view(i) + org_y.h_view(i);
    EXPECT_NEAR_KK(expected, y.h_view(i), 2 * max_error);
  }
}

template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_axpy_mv(int N, int K) {
  using ScalarA    = typename ViewTypeA::value_type;
  using ScalarB    = typename ViewTypeB::value_type;
  using MagnitudeB = typename Kokkos::ArithTraits<ScalarB>::mag_type;

  view_stride_adapter<ViewTypeA> x("X", N, K);
  view_stride_adapter<ViewTypeB> y("Y", N, K);
  view_stride_adapter<ViewTypeB> org_y("Org_Y", N, K);

  ScalarA a                = 3;
  const MagnitudeB eps     = Kokkos::ArithTraits<ScalarB>::epsilon();
  const MagnitudeB max_val = 10;
  const MagnitudeB max_error =
      (static_cast<MagnitudeB>(Kokkos::ArithTraits<ScalarA>::abs(a)) * max_val + max_val) * eps;

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(max_val, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(max_val, randStart, randEnd);
    Kokkos::fill_random(y.d_view, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(org_y.h_base, y.d_base);
  Kokkos::deep_copy(x.h_base, x.d_base);

  KokkosBlas::axpy(a, x.d_view, y.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(static_cast<ScalarB>(a * x.h_view(i, j) + org_y.h_view(i, j)), y.h_view(i, j), 2 * max_error);
    }
  }

  // reset y to orig, and run again with const-valued x
  Kokkos::deep_copy(y.d_base, org_y.h_base);
  KokkosBlas::axpy(a, x.d_view, y.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(static_cast<ScalarB>(a * x.h_view(i, j) + org_y.h_view(i, j)), y.h_view(i, j), 2 * max_error);
    }
  }
}
}  // namespace Test

template <class ScalarA, class ScalarB, class Device>
int test_axpy() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_axpy<view_type_a_ll, view_type_b_ll, Device>(0);
  Test::impl_test_axpy<view_type_a_ll, view_type_b_ll, Device>(13);
  Test::impl_test_axpy<view_type_a_ll, view_type_b_ll, Device>(1024);
  // Test::impl_test_axpy<view_type_a_ll, view_type_b_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_axpy<view_type_a_lr, view_type_b_lr, Device>(0);
  Test::impl_test_axpy<view_type_a_lr, view_type_b_lr, Device>(13);
  Test::impl_test_axpy<view_type_a_lr, view_type_b_lr, Device>(1024);
  // Test::impl_test_axpy<view_type_a_lr, view_type_b_lr, Device>(132231);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_axpy<view_type_a_ls, view_type_b_ls, Device>(0);
  Test::impl_test_axpy<view_type_a_ls, view_type_b_ls, Device>(13);
  Test::impl_test_axpy<view_type_a_ls, view_type_b_ls, Device>(1024);
  // Test::impl_test_axpy<view_type_a_ls, view_type_b_ls, Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_axpy<view_type_a_ls, view_type_b_ll, Device>(1024);
  Test::impl_test_axpy<view_type_a_ll, view_type_b_ls, Device>(1024);
#endif

  return 1;
}

template <class ScalarA, class ScalarB, class Device>
int test_axpy_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_axpy_mv<view_type_a_ll, view_type_b_ll, Device>(0, 5);
  Test::impl_test_axpy_mv<view_type_a_ll, view_type_b_ll, Device>(13, 5);
  Test::impl_test_axpy_mv<view_type_a_ll, view_type_b_ll, Device>(1024, 5);
  // Test::impl_test_axpy_mv<view_type_a_ll, view_type_b_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_axpy_mv<view_type_a_lr, view_type_b_lr, Device>(0, 5);
  Test::impl_test_axpy_mv<view_type_a_lr, view_type_b_lr, Device>(13, 5);
  Test::impl_test_axpy_mv<view_type_a_lr, view_type_b_lr, Device>(1024, 5);
  // Test::impl_test_axpy_mv<view_type_a_lr, view_type_b_lr, Device>(132231,5);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_axpy_mv<view_type_a_ls, view_type_b_ls, Device>(0, 5);
  Test::impl_test_axpy_mv<view_type_a_ls, view_type_b_ls, Device>(13, 5);
  Test::impl_test_axpy_mv<view_type_a_ls, view_type_b_ls, Device>(1024, 5);
  // Test::impl_test_axpy_mv<view_type_a_ls, view_type_b_ls, Device>(132231,5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_axpy_mv<view_type_a_ls, view_type_b_ll, Device>(1024, 5);
  Test::impl_test_axpy_mv<view_type_a_ll, view_type_b_ls, Device>(1024, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpy_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpy_float");
  test_axpy<float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpy_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpy_mv_float");
  test_axpy_mv<float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpy_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpy_double");
  test_axpy<double, double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpy_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpy_mv_double");
  test_axpy_mv<double, double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpy_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpy_complex_double");
  test_axpy<Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpy_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpy_mv_complex_double");
  test_axpy_mv<Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpy_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpy_int");
  test_axpy<int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpy_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpy_mv_int");
  test_axpy_mv<int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, axpy_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpy_double_int");
  test_axpy<double, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpy_double_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpy_mv_double_int");
  test_axpy_mv<double, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
