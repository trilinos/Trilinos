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
#include <KokkosBlas1_iamax.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class Device>
void impl_test_iamax(int N) {
  typedef typename ViewTypeA::non_const_value_type ScalarA;
  typedef Kokkos::ArithTraits<ScalarA> AT;
  typedef typename AT::mag_type mag_type;
  using size_type = typename ViewTypeA::size_type;

  view_stride_adapter<ViewTypeA> a("X", N);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  ScalarA randStart, randEnd;
  Test::getRandomBounds(10.0, randStart, randEnd);
  Kokkos::fill_random(a.d_view, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(a.h_base, a.d_base);

  mag_type expected_result   = Kokkos::ArithTraits<mag_type>::min();
  size_type expected_max_loc = 0;
  for (int i = 0; i < N; i++) {
    mag_type val = AT::abs(a.h_view(i));
    if (val > expected_result) {
      expected_result  = val;
      expected_max_loc = i + 1;
    }
  }

  if (N == 0) {
    expected_result  = typename AT::mag_type(0);
    expected_max_loc = 0;
  }

  {
    // printf("impl_test_iamax -- return result as a scalar on host -- N %d\n",
    // N);
    size_type nonconst_max_loc = KokkosBlas::iamax(a.d_view);
    ASSERT_EQ(nonconst_max_loc, expected_max_loc);

    size_type const_max_loc = KokkosBlas::iamax(a.d_view_const);
    ASSERT_EQ(const_max_loc, expected_max_loc);
  }

  {
    // printf("impl_test_iamax -- return result as a 0-D View on host -- N
    // %d\n", N);
    typedef Kokkos::View<size_type, typename ViewTypeA::array_layout, Kokkos::HostSpace> ViewType0D;
    ViewType0D r("Iamax::Result 0-D View on host", typename ViewTypeA::array_layout());

    KokkosBlas::iamax(r, a.d_view);
    Kokkos::fence();
    size_type nonconst_max_loc = r();
    ASSERT_EQ(nonconst_max_loc, expected_max_loc);

    KokkosBlas::iamax(r, a.d_view_const);
    size_type const_max_loc = r();
    ASSERT_EQ(const_max_loc, expected_max_loc);
  }

  {
    // printf("impl_test_iamax -- return result as a 0-D View on device -- N
    // %d\n", N);
    typedef Kokkos::View<size_type, typename ViewTypeA::array_layout, Device> ViewType0D;
    ViewType0D r("Iamax::Result 0-D View on device", typename ViewTypeA::array_layout());
    typename ViewType0D::HostMirror h_r = Kokkos::create_mirror_view(r);

    size_type nonconst_max_loc, const_max_loc;

    KokkosBlas::iamax(r, a.d_view);
    Kokkos::deep_copy(h_r, r);

    nonconst_max_loc = h_r();

    ASSERT_EQ(nonconst_max_loc, expected_max_loc);

    KokkosBlas::iamax(r, a.d_view_const);
    Kokkos::deep_copy(h_r, r);

    const_max_loc = h_r();

    ASSERT_EQ(const_max_loc, expected_max_loc);
  }
}

template <class ViewTypeA, class Device>
void impl_test_iamax_mv(int N, int K) {
  typedef typename ViewTypeA::non_const_value_type ScalarA;
  typedef Kokkos::ArithTraits<ScalarA> AT;
  typedef typename AT::mag_type mag_type;
  typedef typename ViewTypeA::size_type size_type;

  view_stride_adapter<ViewTypeA> a("A", N, K);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  ScalarA randStart, randEnd;
  Test::getRandomBounds(10.0, randStart, randEnd);
  Kokkos::fill_random(a.d_view, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(a.h_base, a.d_base);

  mag_type* expected_result   = new mag_type[K];
  size_type* expected_max_loc = new size_type[K];

  for (int j = 0; j < K; j++) {
    expected_result[j] = Kokkos::ArithTraits<mag_type>::min();
    for (int i = 0; i < N; i++) {
      mag_type val = AT::abs(a.h_view(i, j));
      if (val > expected_result[j]) {
        expected_result[j]  = val;
        expected_max_loc[j] = i + 1;
      }
    }
    if (N == 0) {
      expected_result[j]  = mag_type(0);
      expected_max_loc[j] = size_type(0);
    }
  }

  {
    // printf("impl_test_iamax_mv -- return results as a 1-D View on host -- N
    // %d\n", N);
    Kokkos::View<size_type*, Kokkos::HostSpace> rcontig("Iamax::Result View on host", K);
    Kokkos::View<size_type*, typename ViewTypeA::array_layout, Kokkos::HostSpace> r = rcontig;

    KokkosBlas::iamax(r, a.d_view);
    Kokkos::fence();

    for (int k = 0; k < K; k++) {
      size_type nonconst_result = r(k);
      size_type exp_result      = expected_max_loc[k];
      ASSERT_EQ(nonconst_result, exp_result);
    }

    KokkosBlas::iamax(r, a.d_view_const);
    Kokkos::fence();

    for (int k = 0; k < K; k++) {
      size_type const_result = r(k);
      size_type exp_result   = expected_max_loc[k];
      ASSERT_EQ(const_result, exp_result);
    }
  }

  {
    // printf("impl_test_iamax_mv -- return results as a 1-D View on device -- N
    // %d\n", N);
    Kokkos::View<size_type*, Device> rcontig("Iamax::Result View on host", K);
    Kokkos::View<size_type*, typename ViewTypeA::array_layout, Device> r = rcontig;
    typename Kokkos::View<size_type*, typename ViewTypeA::array_layout, Device>::HostMirror h_r =
        Kokkos::create_mirror_view(rcontig);

    KokkosBlas::iamax(r, a.d_view);
    Kokkos::deep_copy(h_r, r);

    for (int k = 0; k < K; k++) {
      size_type nonconst_result = h_r(k);
      size_type exp_result      = expected_max_loc[k];
      ASSERT_EQ(nonconst_result, exp_result);
    }

    KokkosBlas::iamax(r, a.d_view_const);
    Kokkos::deep_copy(h_r, r);

    for (int k = 0; k < K; k++) {
      size_type const_result = h_r(k);
      size_type exp_result   = expected_max_loc[k];
      ASSERT_EQ(const_result, exp_result);
    }
  }

  delete[] expected_result;
  delete[] expected_max_loc;
}
}  // namespace Test

template <class ScalarA, class Device>
int test_iamax() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_iamax<view_type_a_ll, Device>(0);
  Test::impl_test_iamax<view_type_a_ll, Device>(13);
  Test::impl_test_iamax<view_type_a_ll, Device>(1024);
  // Test::impl_test_iamax<view_type_a_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_iamax<view_type_a_lr, Device>(0);
  Test::impl_test_iamax<view_type_a_lr, Device>(13);
  Test::impl_test_iamax<view_type_a_lr, Device>(1024);
  // Test::impl_test_iamax<view_type_a_lr, Device>(132231);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_iamax<view_type_a_ls, Device>(0);
  Test::impl_test_iamax<view_type_a_ls, Device>(13);
  Test::impl_test_iamax<view_type_a_ls, Device>(1024);
  // Test::impl_test_iamax<view_type_a_ls, Device>(132231);
#endif

  return 1;
}

template <class ScalarA, class Device>
int test_iamax_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_iamax_mv<view_type_a_ll, Device>(0, 5);
  Test::impl_test_iamax_mv<view_type_a_ll, Device>(13, 5);
  Test::impl_test_iamax_mv<view_type_a_ll, Device>(1024, 5);
  // Test::impl_test_iamax_mv<view_type_a_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_iamax_mv<view_type_a_lr, Device>(0, 5);
  Test::impl_test_iamax_mv<view_type_a_lr, Device>(13, 5);
  Test::impl_test_iamax_mv<view_type_a_lr, Device>(1024, 5);
  // Test::impl_test_iamax_mv<view_type_a_lr, Device>(132231,5);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_iamax_mv<view_type_a_ls, Device>(0, 5);
  Test::impl_test_iamax_mv<view_type_a_ls, Device>(13, 5);
  Test::impl_test_iamax_mv<view_type_a_ls, Device>(1024, 5);
  // Test::impl_test_iamax_mv<view_type_a_ls, Device>(132231,5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, iamax_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::iamax_float");
  test_iamax<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, iamax_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::iamax_mvfloat");
  test_iamax_mv<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, iamax_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::iamax_double");
  test_iamax<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, iamax_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::iamax_mv_double");
  test_iamax_mv<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, iamax_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::iamax_complex_double");
  test_iamax<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, iamax_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::iamax_mv_complex_double");
  test_iamax_mv<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, iamax_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::iamax_int");
  test_iamax<int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, iamax_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::iamax_mv_int");
  test_iamax_mv<int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
