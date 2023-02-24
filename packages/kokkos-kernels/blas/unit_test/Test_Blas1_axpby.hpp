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
void impl_test_axpby(int N) {
  using ScalarA    = typename ViewTypeA::value_type;
  using ScalarB    = typename ViewTypeB::value_type;
  using MagnitudeB = typename Kokkos::ArithTraits<ScalarB>::mag_type;

  using BaseTypeA = Kokkos::View<
      ScalarA * [2],
      typename std::conditional<std::is_same<typename ViewTypeA::array_layout,
                                             Kokkos::LayoutStride>::value,
                                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,
      Device>;
  using BaseTypeB = Kokkos::View<
      ScalarB * [2],
      typename std::conditional<std::is_same<typename ViewTypeB::array_layout,
                                             Kokkos::LayoutStride>::value,
                                Kokkos::LayoutRight, Kokkos::LayoutLeft>::type,
      Device>;

  ScalarA a = 3;
  ScalarB b = 5;
  // eps should probably be based on ScalarB since that is the type
  // in which the result is computed.
  const MagnitudeB eps     = Kokkos::ArithTraits<ScalarB>::epsilon();
  const MagnitudeB max_val = 10;
  const MagnitudeB max_error =
      (static_cast<MagnitudeB>(Kokkos::ArithTraits<ScalarA>::abs(a)) +
       Kokkos::ArithTraits<ScalarB>::abs(b)) *
      max_val * eps;

  BaseTypeA b_x("X", N);
  BaseTypeB b_y("Y", N);
  BaseTypeB b_org_y("Org_Y", N);

  auto h_b_org_y =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b_org_y);
  ViewTypeA x                        = Kokkos::subview(b_x, Kokkos::ALL(), 0);
  ViewTypeB y                        = Kokkos::subview(b_y, Kokkos::ALL(), 0);
  typename ViewTypeA::const_type c_x = x;
  typename ViewTypeB::const_type c_y = y;

  typename BaseTypeA::HostMirror h_b_x = Kokkos::create_mirror_view(b_x);
  typename BaseTypeB::HostMirror h_b_y = Kokkos::create_mirror_view(b_y);

  typename ViewTypeA::HostMirror h_x = Kokkos::subview(h_b_x, Kokkos::ALL(), 0);
  typename ViewTypeB::HostMirror h_y = Kokkos::subview(h_b_y, Kokkos::ALL(), 0);

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(max_val, randStart, randEnd);
    Kokkos::fill_random(b_x, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(max_val, randStart, randEnd);
    Kokkos::fill_random(b_y, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(b_org_y, b_y);
  Kokkos::deep_copy(h_b_org_y, b_org_y);

  Kokkos::deep_copy(h_b_x, b_x);

  // Run with non-const input (x) and verify
  KokkosBlas::axpby(a, x, b, y);
  Kokkos::deep_copy(h_b_y, b_y);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(static_cast<ScalarB>(a * h_x(i) + b * h_b_org_y(i, 0)),
                   h_y(i), 2 * max_error);
  }

  Kokkos::deep_copy(b_y, b_org_y);
  // Run again with const input (c_x)
  KokkosBlas::axpby(a, c_x, b, y);
  Kokkos::deep_copy(h_b_y, b_y);
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR_KK(static_cast<ScalarB>(a * h_x(i) + b * h_b_org_y(i, 0)),
                   h_y(i), 2 * max_error);
  }
}

template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_axpby_mv(int N, int K) {
  using ScalarA    = typename ViewTypeA::value_type;
  using ScalarB    = typename ViewTypeB::value_type;
  using MagnitudeB = typename Kokkos::ArithTraits<ScalarB>::mag_type;

  typedef multivector_layout_adapter<ViewTypeA> vfA_type;
  typedef multivector_layout_adapter<ViewTypeB> vfB_type;

  typename vfA_type::BaseType b_x("A", N, K);
  typename vfB_type::BaseType b_y("B", N, K);
  typename vfB_type::BaseType b_org_y("B", N, K);

  ViewTypeA x = vfA_type::view(b_x);
  ViewTypeB y = vfB_type::view(b_y);

  typedef multivector_layout_adapter<typename ViewTypeA::HostMirror> h_vfA_type;
  typedef multivector_layout_adapter<typename ViewTypeB::HostMirror> h_vfB_type;

  typename h_vfA_type::BaseType h_b_x = Kokkos::create_mirror_view(b_x);
  typename h_vfB_type::BaseType h_b_y = Kokkos::create_mirror_view(b_y);

  typename ViewTypeA::HostMirror h_x = h_vfA_type::view(h_b_x);
  typename ViewTypeB::HostMirror h_y = h_vfB_type::view(h_b_y);

  ScalarA a                = 3;
  ScalarB b                = 5;
  const MagnitudeB eps     = Kokkos::ArithTraits<ScalarB>::epsilon();
  const MagnitudeB max_val = 10;
  const MagnitudeB max_error =
      (static_cast<MagnitudeB>(Kokkos::ArithTraits<ScalarA>::abs(a)) +
       Kokkos::ArithTraits<ScalarB>::abs(b)) *
      max_val * eps;

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(
      13718);

  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(b_x, rand_pool, randStart, randEnd);
  }
  {
    ScalarB randStart, randEnd;
    Test::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::fill_random(b_y, rand_pool, randStart, randEnd);
  }

  Kokkos::deep_copy(b_org_y, b_y);
  ViewTypeB org_y = vfB_type::view(b_org_y);
  auto h_org_y =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), org_y);

  Kokkos::deep_copy(h_b_x, b_x);
  Kokkos::deep_copy(h_b_y, b_y);

  typename ViewTypeA::const_type c_x = x;

  Kokkos::View<ScalarB*, Kokkos::HostSpace> r("Dot::Result", K);

  KokkosBlas::axpby(a, x, b, y);
  Kokkos::deep_copy(h_b_y, b_y);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(static_cast<ScalarB>(a * h_x(i, j) + b * h_org_y(i, j)),
                     h_y(i, j), 2 * max_error);
    }
  }

  Kokkos::deep_copy(b_y, b_org_y);
  KokkosBlas::axpby(a, c_x, b, y);
  Kokkos::deep_copy(h_b_y, b_y);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      EXPECT_NEAR_KK(static_cast<ScalarB>(a * h_x(i, j) + b * h_org_y(i, j)),
                     h_y(i, j), 2 * max_error);
    }
  }
}
}  // namespace Test

template <class ScalarA, class ScalarB, class Device>
int test_axpby() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_axpby<view_type_a_ll, view_type_b_ll, Device>(0);
  Test::impl_test_axpby<view_type_a_ll, view_type_b_ll, Device>(13);
  Test::impl_test_axpby<view_type_a_ll, view_type_b_ll, Device>(1024);
  Test::impl_test_axpby<view_type_a_ll, view_type_b_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_axpby<view_type_a_lr, view_type_b_lr, Device>(0);
  Test::impl_test_axpby<view_type_a_lr, view_type_b_lr, Device>(13);
  Test::impl_test_axpby<view_type_a_lr, view_type_b_lr, Device>(1024);
  Test::impl_test_axpby<view_type_a_lr, view_type_b_lr, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB*, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_axpby<view_type_a_ls, view_type_b_ls, Device>(0);
  Test::impl_test_axpby<view_type_a_ls, view_type_b_ls, Device>(13);
  Test::impl_test_axpby<view_type_a_ls, view_type_b_ls, Device>(1024);
  Test::impl_test_axpby<view_type_a_ls, view_type_b_ls, Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_axpby<view_type_a_ls, view_type_b_ll, Device>(1024);
  Test::impl_test_axpby<view_type_a_ll, view_type_b_ls, Device>(1024);
#endif

  return 1;
}

template <class ScalarA, class ScalarB, class Device>
int test_axpby_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&      \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_axpby_mv<view_type_a_ll, view_type_b_ll, Device>(0, 5);
  Test::impl_test_axpby_mv<view_type_a_ll, view_type_b_ll, Device>(13, 5);
  Test::impl_test_axpby_mv<view_type_a_ll, view_type_b_ll, Device>(1024, 5);
  Test::impl_test_axpby_mv<view_type_a_ll, view_type_b_ll, Device>(132231, 5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&       \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_axpby_mv<view_type_a_lr, view_type_b_lr, Device>(0, 5);
  Test::impl_test_axpby_mv<view_type_a_lr, view_type_b_lr, Device>(13, 5);
  Test::impl_test_axpby_mv<view_type_a_lr, view_type_b_lr, Device>(1024, 5);
  Test::impl_test_axpby_mv<view_type_a_lr, view_type_b_lr, Device>(132231, 5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&        \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB**, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_axpby_mv<view_type_a_ls, view_type_b_ls, Device>(0, 5);
  Test::impl_test_axpby_mv<view_type_a_ls, view_type_b_ls, Device>(13, 5);
  Test::impl_test_axpby_mv<view_type_a_ls, view_type_b_ls, Device>(1024, 5);
  Test::impl_test_axpby_mv<view_type_a_ls, view_type_b_ls, Device>(132231, 5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_axpby_mv<view_type_a_ls, view_type_b_ll, Device>(1024, 5);
  Test::impl_test_axpby_mv<view_type_a_ll, view_type_b_ls, Device>(1024, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpby_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_float");
  test_axpby<float, float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpby_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_mv_float");
  test_axpby_mv<float, float, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&  \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpby_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_double");
  test_axpby<double, double, TestExecSpace>();
}
TEST_F(TestCategory, axpby_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_mv_double");
  test_axpby_mv<double, double, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) &&          \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpby_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_complex_double");
  test_axpby<Kokkos::complex<double>, Kokkos::complex<double>, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpby_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_mv_complex_double");
  test_axpby_mv<Kokkos::complex<double>, Kokkos::complex<double>,
                TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) ||   \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, axpby_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_int");
  test_axpby<int, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpby_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_mv_int");
  test_axpby_mv<int, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && \
    !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, axpby_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_double_int");
  test_axpby<double, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, axpby_double_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::axpby_mv_double_int");
  test_axpby_mv<double, int, TestExecSpace>();
  Kokkos::Profiling::popRegion();
}
#endif
