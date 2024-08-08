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
#include <KokkosBlas1_scal.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_scal(int N) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef Kokkos::ArithTraits<ScalarA> AT;

  ScalarA a(3);
  typename AT::mag_type eps = AT::epsilon() * 1000;

  view_stride_adapter<ViewTypeA> x("X", N);
  view_stride_adapter<ViewTypeB> y("Y", N);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(x.h_base, x.d_base);

  KokkosBlas::scal(y.d_view, a, x.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(static_cast<ScalarB>(a * x.h_view(i)), y.h_view(i), eps);
  }

  // Zero out y again and run with const input
  Kokkos::deep_copy(y.d_view, Kokkos::ArithTraits<ScalarB>::zero());
  KokkosBlas::scal(y.d_view, a, x.d_view_const);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(static_cast<ScalarB>(a * x.h_view(i)), y.h_view(i), eps);
  }
}

template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_scal_mv(int N, int K) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef Kokkos::ArithTraits<ScalarA> AT;

  view_stride_adapter<ViewTypeA> x("X", N, K);
  view_stride_adapter<ViewTypeB> y("Y", N, K);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(1.0, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(x.h_base, x.d_base);

  ScalarA a(3.0);

  typename AT::mag_type eps = AT::epsilon() * 1000;

  Kokkos::View<ScalarB*, Kokkos::HostSpace> r("Dot::Result", K);

  KokkosBlas::scal(y.d_view, a, x.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(static_cast<ScalarB>(a * x.h_view(i, j)), y.h_view(i, j), eps);
    }
  }

  // Zero out y again, and run again with const input
  Kokkos::deep_copy(y.d_view, Kokkos::ArithTraits<ScalarB>::zero());
  KokkosBlas::scal(y.d_view, a, x.d_view_const);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(static_cast<ScalarB>(a * x.h_view(i, j)), y.h_view(i, j), eps);
    }
  }

  // Generate 'params' view with dimension == number of multivectors; each entry
  // will be different scalar to scale y
  Kokkos::View<ScalarA*, Device> params("Params", K);
  for (int j = 0; j < K; j++) {
    Kokkos::View<ScalarA, Device> param_j(params, j);
    Kokkos::deep_copy(param_j, ScalarA(3 + j));
  }

  auto h_params = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), params);

  Kokkos::deep_copy(y.d_view, Kokkos::ArithTraits<ScalarB>::zero());
  KokkosBlas::scal(y.d_view, params, x.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(static_cast<ScalarB>(h_params(j) * x.h_view(i, j)), y.h_view(i, j), eps);
    }
  }

  Kokkos::deep_copy(y.d_view, Kokkos::ArithTraits<ScalarB>::zero());
  KokkosBlas::scal(y.d_view, params, x.d_view_const);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(static_cast<ScalarB>(h_params(j) * x.h_view(i, j)), y.h_view(i, j), eps);
    }
  }
}
}  // namespace Test

template <class ScalarA, class ScalarB, class Device>
int test_scal() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_scal<view_type_a_ll, view_type_b_ll, Device>(0);
  Test::impl_test_scal<view_type_a_ll, view_type_b_ll, Device>(13);
  Test::impl_test_scal<view_type_a_ll, view_type_b_ll, Device>(1024);
  // Test::impl_test_scal<view_type_a_ll, view_type_b_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_scal<view_type_a_lr, view_type_b_lr, Device>(0);
  Test::impl_test_scal<view_type_a_lr, view_type_b_lr, Device>(13);
  Test::impl_test_scal<view_type_a_lr, view_type_b_lr, Device>(1024);
  // Test::impl_test_scal<view_type_a_lr, view_type_b_lr, Device>(132231);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_scal<view_type_a_ls, view_type_b_ls, Device>(0);
  Test::impl_test_scal<view_type_a_ls, view_type_b_ls, Device>(13);
  Test::impl_test_scal<view_type_a_ls, view_type_b_ls, Device>(1024);
  // Test::impl_test_scal<view_type_a_ls, view_type_b_ls, Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_scal<view_type_a_ls, view_type_b_ll, Device>(1024);
  Test::impl_test_scal<view_type_a_ll, view_type_b_ls, Device>(1024);
#endif

  return 1;
}

template <class ScalarA, class ScalarB, class Device>
int test_scal_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_scal_mv<view_type_a_ll, view_type_b_ll, Device>(0, 5);
  Test::impl_test_scal_mv<view_type_a_ll, view_type_b_ll, Device>(13, 5);
  Test::impl_test_scal_mv<view_type_a_ll, view_type_b_ll, Device>(1024, 5);
  // Test::impl_test_scal_mv<view_type_a_ll, view_type_b_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_scal_mv<view_type_a_lr, view_type_b_lr, Device>(0, 5);
  Test::impl_test_scal_mv<view_type_a_lr, view_type_b_lr, Device>(13, 5);
  Test::impl_test_scal_mv<view_type_a_lr, view_type_b_lr, Device>(1024, 5);
  // Test::impl_test_scal_mv<view_type_a_lr, view_type_b_lr, Device>(132231,5);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_scal_mv<view_type_a_ls, view_type_b_ls, Device>(0, 5);
  Test::impl_test_scal_mv<view_type_a_ls, view_type_b_ls, Device>(13, 5);
  Test::impl_test_scal_mv<view_type_a_ls, view_type_b_ls, Device>(1024, 5);
  // Test::impl_test_scal_mv<view_type_a_ls, view_type_b_ls, Device>(132231,5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_scal_mv<view_type_a_ls, view_type_b_ll, Device>(1024, 5);
  Test::impl_test_scal_mv<view_type_a_ll, view_type_b_ls, Device>(1024, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, scal_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_float");
  test_scal<float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, scal_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_mv_float");
  test_scal_mv<float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, scal_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_double");
  test_scal<double, double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, scal_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_mv_double");
  test_scal_mv<double, double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, scal_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_complex_double");
  test_scal<Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, scal_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_mv_complex_double");
  test_scal_mv<Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, scal_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_int");
  test_scal<int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, scal_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_mv_int");
  test_scal_mv<int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, scal_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_double_int");
  test_scal<double, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, scal_mv_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::scal_mv_double_int");
  test_scal_mv<double, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
