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
#include <KokkosBlas1_reciprocal.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_reciprocal(int N) {
  using ScalarA    = typename ViewTypeA::value_type;
  using ScalarB    = typename ViewTypeB::value_type;
  using AT         = Kokkos::ArithTraits<ScalarA>;
  using MagnitudeA = typename AT::mag_type;
  using MagnitudeB = typename Kokkos::ArithTraits<ScalarB>::mag_type;

  const MagnitudeB eps     = Kokkos::ArithTraits<ScalarB>::epsilon();
  const MagnitudeA one     = AT::abs(AT::one());
  const MagnitudeA max_val = 10;

  view_stride_adapter<ViewTypeA> x("X", N);
  view_stride_adapter<ViewTypeB> y("Y", N);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(max_val, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, one, randEnd);
  }

  Kokkos::deep_copy(x.h_base, x.d_base);

  KokkosBlas::reciprocal(y.d_view, x.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; ++i) {
    EXPECT_NEAR_KK(y.h_view(i), ScalarB(one / x.h_view(i)), 2 * eps);
  }

  // Zero out y again, and run again with const input
  Kokkos::deep_copy(y.d_view, Kokkos::ArithTraits<ScalarB>::zero());

  KokkosBlas::reciprocal(y.d_view, x.d_view_const);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int i = 0; i < N; ++i) {
    EXPECT_NEAR_KK(y.h_view(i), ScalarB(one / x.h_view(i)), 2 * eps);
  }
}

template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_reciprocal_mv(int N, int K) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;

  view_stride_adapter<ViewTypeA> x("X", N, K);
  view_stride_adapter<ViewTypeB> y("Y", N, K);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(10, randStart, randEnd);
    Kokkos::fill_random(x.d_view, rand_pool, Kokkos::ArithTraits<ScalarA>::one(), randEnd);
  }

  Kokkos::deep_copy(x.h_base, x.d_base);

  KokkosBlas::reciprocal(y.d_view, x.d_view);

  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int j = 0; j < K; ++j) {
    for (int i = 0; i < N; ++i) {
      EXPECT_NEAR_KK(y.h_view(i, j), Kokkos::ArithTraits<ScalarB>::one() / ScalarB(x.h_view(i, j)),
                     2 * Kokkos::ArithTraits<ScalarB>::epsilon());
    }
  }

  // Zero out y again, and run again with const input
  Kokkos::deep_copy(y.d_view, Kokkos::ArithTraits<ScalarB>::zero());

  KokkosBlas::reciprocal(y.d_view, x.d_view_const);
  Kokkos::deep_copy(y.h_base, y.d_base);
  for (int j = 0; j < K; j++) {
    for (int i = 0; i < N; ++i) {
      EXPECT_NEAR_KK(y.h_view(i, j), Kokkos::ArithTraits<ScalarB>::one() / ScalarB(x.h_view(i, j)),
                     2 * Kokkos::ArithTraits<ScalarB>::epsilon());
    }
  }
}
}  // namespace Test

template <class ScalarA, class ScalarB, class Device>
int test_reciprocal() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_reciprocal<view_type_a_ll, view_type_b_ll, Device>(0);
  Test::impl_test_reciprocal<view_type_a_ll, view_type_b_ll, Device>(13);
  Test::impl_test_reciprocal<view_type_a_ll, view_type_b_ll, Device>(1024);
  // Test::impl_test_reciprocal<view_type_a_ll, view_type_b_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_reciprocal<view_type_a_lr, view_type_b_lr, Device>(0);
  Test::impl_test_reciprocal<view_type_a_lr, view_type_b_lr, Device>(13);
  Test::impl_test_reciprocal<view_type_a_lr, view_type_b_lr, Device>(1024);
  // Test::impl_test_reciprocal<view_type_a_lr, view_type_b_lr, Device>(132231);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_reciprocal<view_type_a_ls, view_type_b_ls, Device>(0);
  Test::impl_test_reciprocal<view_type_a_ls, view_type_b_ls, Device>(13);
  Test::impl_test_reciprocal<view_type_a_ls, view_type_b_ls, Device>(1024);
  // Test::impl_test_reciprocal<view_type_a_ls, view_type_b_ls, Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_reciprocal<view_type_a_ls, view_type_b_ll, Device>(1024);
  Test::impl_test_reciprocal<view_type_a_ll, view_type_b_ls, Device>(1024);
#endif

  return 1;
}

template <class ScalarA, class ScalarB, class Device>
int test_reciprocal_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_reciprocal_mv<view_type_a_ll, view_type_b_ll, Device>(0, 5);
  Test::impl_test_reciprocal_mv<view_type_a_ll, view_type_b_ll, Device>(13, 5);
  Test::impl_test_reciprocal_mv<view_type_a_ll, view_type_b_ll, Device>(1024, 5);
  // Test::impl_test_reciprocal_mv<view_type_a_ll, view_type_b_ll,
  // Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_reciprocal_mv<view_type_a_lr, view_type_b_lr, Device>(0, 5);
  Test::impl_test_reciprocal_mv<view_type_a_lr, view_type_b_lr, Device>(13, 5);
  Test::impl_test_reciprocal_mv<view_type_a_lr, view_type_b_lr, Device>(1024, 5);
  // Test::impl_test_reciprocal_mv<view_type_a_lr, view_type_b_lr,
  // Device>(132231,5);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_reciprocal_mv<view_type_a_ls, view_type_b_ls, Device>(0, 5);
  Test::impl_test_reciprocal_mv<view_type_a_ls, view_type_b_ls, Device>(13, 5);
  Test::impl_test_reciprocal_mv<view_type_a_ls, view_type_b_ls, Device>(1024, 5);
  // Test::impl_test_reciprocal_mv<view_type_a_ls, view_type_b_ls,
  // Device>(132231,5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_reciprocal_mv<view_type_a_ls, view_type_b_ll, Device>(1024, 5);
  Test::impl_test_reciprocal_mv<view_type_a_ll, view_type_b_ls, Device>(1024, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, reciprocal_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::reciprocal_float");
  test_reciprocal<float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, reciprocal_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::reciprocal_mv_float");
  test_reciprocal_mv<float, float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, reciprocal_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::reciprocal_double");
  test_reciprocal<double, double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, reciprocal_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::reciprocal_mv_double");
  test_reciprocal_mv<double, double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, reciprocal_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::reciprocal_complex_double");
  test_reciprocal<Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, reciprocal_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::reciprocal_mv_complex_double");
  test_reciprocal_mv<Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, reciprocal_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::reciprocal_int");
  test_reciprocal<int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, reciprocal_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::reciprocal_mv_int");
  test_reciprocal_mv<int, int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

/*
#if !defined(KOKKOSKERNELS_ETI_ONLY) &&
!defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS) TEST_F( TestCategory,
reciprocal_double_int ) { test_reciprocal<double,int,TestDevice> ();
}
TEST_F( TestCategory, reciprocal_double_mv_int ) {
    test_reciprocal_mv<double,int,TestDevice> ();
}
#endif
*/
