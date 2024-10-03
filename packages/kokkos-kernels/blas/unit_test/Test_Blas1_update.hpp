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
#include <KokkosBlas1_update.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class Device>
void impl_test_update(int N) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef typename ViewTypeC::value_type ScalarC;

  ScalarA a  = 3;
  ScalarB b  = 5;
  ScalarC c  = 7;
  double eps = std::is_same<ScalarC, float>::value ? 2 * 1e-5 : 1e-7;

  view_stride_adapter<ViewTypeA> x("X", N);
  view_stride_adapter<ViewTypeB> y("Y", N);
  view_stride_adapter<ViewTypeC> z("Z", N);
  view_stride_adapter<ViewTypeC> org_z("Org_Z", N);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(y.d_view, rand_pool, randStart, randEnd);
  }
  {
    ScalarC randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(z.d_view, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(org_z.h_base, z.d_base);

  Kokkos::deep_copy(x.h_base, x.d_base);
  Kokkos::deep_copy(y.h_base, y.d_base);

  KokkosBlas::update(a, x.d_view, b, y.d_view, c, z.d_view);
  Kokkos::deep_copy(z.h_base, z.d_base);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(static_cast<ScalarC>(a * x.h_view(i) + b * y.h_view(i) + c * org_z.h_view(i)), z.h_view(i), eps);
  }

  Kokkos::deep_copy(z.d_base, org_z.h_base);
  KokkosBlas::update(a, x.d_view_const, b, y.d_view, c, z.d_view);
  Kokkos::deep_copy(z.h_base, z.d_base);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(static_cast<ScalarC>(a * x.h_view(i) + b * y.h_view(i) + c * org_z.h_view(i)), z.h_view(i), eps);
  }

  Kokkos::deep_copy(z.d_base, org_z.h_base);
  KokkosBlas::update(a, x.d_view_const, b, y.d_view_const, c, z.d_view);
  Kokkos::deep_copy(z.h_base, z.d_base);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(static_cast<ScalarC>(a * x.h_view(i) + b * y.h_view(i) + c * org_z.h_view(i)), z.h_view(i), eps);
  }
}

template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class Device>
void impl_test_update_mv(int N, int K) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef typename ViewTypeC::value_type ScalarC;

  view_stride_adapter<ViewTypeA> x("X", N, K);
  view_stride_adapter<ViewTypeB> y("Y", N, K);
  view_stride_adapter<ViewTypeC> z("Z", N, K);
  view_stride_adapter<ViewTypeC> org_z("Org_Z", N, K);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(y.d_view, rand_pool, randStart, randEnd);
  }
  {
    ScalarC randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(z.d_view, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(org_z.h_base, z.d_base);

  Kokkos::deep_copy(x.h_base, x.d_base);
  Kokkos::deep_copy(y.h_base, y.d_base);

  ScalarA a = 3;
  ScalarB b = 5;
  ScalarC c = 5;

  double eps = std::is_same<ScalarA, float>::value ? 2 * 1e-5 : 1e-7;

  KokkosBlas::update(a, x.d_view, b, y.d_view, c, z.d_view);
  Kokkos::deep_copy(z.h_base, z.d_base);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(static_cast<ScalarC>(a * x.h_view(i, j) + b * y.h_view(i, j) + c * org_z.h_view(i, j)),
                     z.h_view(i, j), eps);
    }
  }

  Kokkos::deep_copy(z.d_base, org_z.h_base);
  KokkosBlas::update(a, x.d_view_const, b, y.d_view, c, z.d_view);
  Kokkos::deep_copy(z.h_base, z.d_base);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(static_cast<ScalarC>(a * x.h_view(i, j) + b * y.h_view(i, j) + c * org_z.h_view(i, j)),
                     z.h_view(i, j), eps);
    }
  }
}
}  // namespace Test

template <class ScalarA, class ScalarB, class ScalarC, class Device>
int test_update() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  typedef Kokkos::View<ScalarC*, Kokkos::LayoutLeft, Device> view_type_c_ll;
  Test::impl_test_update<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(0);
  Test::impl_test_update<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(13);
  Test::impl_test_update<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(1024);
  // Test::impl_test_update<view_type_a_ll, view_type_b_ll, view_type_c_ll,
  // Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutRight, Device> view_type_b_lr;
  typedef Kokkos::View<ScalarC*, Kokkos::LayoutRight, Device> view_type_c_lr;
  Test::impl_test_update<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(0);
  Test::impl_test_update<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(13);
  Test::impl_test_update<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(1024);
  // Test::impl_test_update<view_type_a_lr, view_type_b_lr, view_type_c_lr,
  // Device>(132231);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, Device> view_type_b_ls;
  typedef Kokkos::View<ScalarC*, Kokkos::LayoutStride, Device> view_type_c_ls;
  Test::impl_test_update<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(0);
  Test::impl_test_update<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(13);
  Test::impl_test_update<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(1024);
  // Test::impl_test_update<view_type_a_ls, view_type_b_ls, view_type_c_ls,
  // Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_update<view_type_a_ls, view_type_b_ll, view_type_c_lr, Device>(1024);
  Test::impl_test_update<view_type_a_ll, view_type_b_ls, view_type_c_lr, Device>(1024);
#endif

  return 1;
}

template <class ScalarA, class ScalarB, class ScalarC, class Device>
int test_update_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  typedef Kokkos::View<ScalarC**, Kokkos::LayoutLeft, Device> view_type_c_ll;
  Test::impl_test_update_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(0, 5);
  Test::impl_test_update_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(13, 5);
  Test::impl_test_update_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(1024, 5);
  Test::impl_test_update_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(132231, 5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  typedef Kokkos::View<ScalarC**, Kokkos::LayoutRight, Device> view_type_c_lr;
  Test::impl_test_update_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(0, 5);
  Test::impl_test_update_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(13, 5);
  Test::impl_test_update_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(1024, 5);
  Test::impl_test_update_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(132231, 5);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutStride, Device> view_type_b_ls;
  typedef Kokkos::View<ScalarC**, Kokkos::LayoutStride, Device> view_type_c_ls;
  Test::impl_test_update_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(0, 5);
  Test::impl_test_update_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(13, 5);
  Test::impl_test_update_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(1024, 5);
  Test::impl_test_update_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(132231, 5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_update_mv<view_type_a_ls, view_type_b_ll, view_type_c_lr, Device>(1024, 5);
  Test::impl_test_update_mv<view_type_a_ll, view_type_b_ls, view_type_c_lr, Device>(1024, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, update_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::update_float");
  test_update<float, float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, update_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::update_mv_float");
  test_update_mv<float, float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, update_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::update_double");
  test_update<double, double, double, TestDevice>();
}
TEST_F(TestCategory, update_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::update_mv_double");
  test_update_mv<double, double, double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, update_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::update_complex_double");
  test_update<Kokkos::complex<double>, Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, update_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::update_mv_complex_double");
  test_update_mv<Kokkos::complex<double>, Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, update_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::update_int");
  test_update<int, int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, update_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::update_mv_int");
  test_update_mv<int, int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, update_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::update_double_int");
  test_update<double, int, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, update_mv_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::update_mv_double_int");
  test_update_mv<double, int, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
