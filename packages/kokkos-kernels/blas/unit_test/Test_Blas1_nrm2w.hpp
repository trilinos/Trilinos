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
#include <KokkosBlas1_nrm2w.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class Device>
void impl_test_nrm2w(int N) {
  using ScalarA    = typename ViewTypeA::value_type;
  using AT         = Kokkos::ArithTraits<ScalarA>;
  using MagnitudeA = typename AT::mag_type;

  view_stride_adapter<ViewTypeA> a("A", N);
  view_stride_adapter<ViewTypeA> w("W", N);

  constexpr MagnitudeA max_val = 10;
  const MagnitudeA eps         = AT::epsilon();
  const MagnitudeA max_error   = max_val * std::sqrt(static_cast<MagnitudeA>(N)) * eps;

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  ScalarA randStart, randEnd;
  Test::getRandomBounds(max_val, randStart, randEnd);
  Kokkos::fill_random(a.d_view, rand_pool, randStart, randEnd);
  Kokkos::fill_random(w.d_view, rand_pool, AT::one(),
                      randEnd);  // Avoid divide by 0

  Kokkos::deep_copy(a.h_base, a.d_base);
  Kokkos::deep_copy(w.h_base, w.d_base);

  typename AT::mag_type expected_result = 0;
  for (int i = 0; i < N; i++) {
    typename AT::mag_type term = AT::abs(a.h_view(i)) / AT::abs(w.h_view(i));
    expected_result += term * term;
  }
  expected_result = Kokkos::ArithTraits<typename AT::mag_type>::sqrt(expected_result);

  typename AT::mag_type nonconst_result = KokkosBlas::nrm2w(a.d_view, w.d_view);
  EXPECT_NEAR_KK(nonconst_result, expected_result, max_error);
}

template <class ViewTypeA, class Device>
void impl_test_nrm2w_mv(int N, int K) {
  using ScalarA    = typename ViewTypeA::value_type;
  using AT         = Kokkos::ArithTraits<ScalarA>;
  using MagnitudeA = typename AT::mag_type;

  view_stride_adapter<ViewTypeA> a("A", N, K);
  view_stride_adapter<ViewTypeA> w("W", N, K);

  constexpr MagnitudeA max_val = 10;
  const MagnitudeA eps         = AT::epsilon();
  const MagnitudeA max_error   = max_val * std::sqrt(static_cast<MagnitudeA>(N)) * eps;

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  ScalarA randStart, randEnd;
  Test::getRandomBounds(max_val, randStart, randEnd);
  Kokkos::fill_random(a.d_view, rand_pool, randStart, randEnd);
  Kokkos::fill_random(w.d_view, rand_pool, AT::one(),
                      randEnd);  // Avoid dividing by 0

  Kokkos::deep_copy(a.h_base, a.d_base);
  Kokkos::deep_copy(w.h_base, w.d_base);

  typename AT::mag_type* expected_result = new typename AT::mag_type[K];
  for (int j = 0; j < K; j++) {
    expected_result[j] = typename AT::mag_type();
    for (int i = 0; i < N; i++) {
      typename AT::mag_type term = AT::abs(a.h_view(i, j)) / AT::abs(w.h_view(i, j));
      expected_result[j] += term * term;
    }
    expected_result[j] = Kokkos::ArithTraits<typename AT::mag_type>::sqrt(expected_result[j]);
  }

  Kokkos::View<typename AT::mag_type*, Device> r("Dot::Result", K);
  KokkosBlas::nrm2w(r, a.d_view, w.d_view);
  auto r_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), r);

  for (int k = 0; k < K; k++) {
    typename AT::mag_type nonconst_result = r_host(k);
    EXPECT_NEAR_KK(nonconst_result, expected_result[k], max_error);
  }

  delete[] expected_result;
}
}  // namespace Test

template <class ScalarA, class Device>
int test_nrm2w() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrm2w<view_type_a_ll, Device>(0);
  Test::impl_test_nrm2w<view_type_a_ll, Device>(13);
  Test::impl_test_nrm2w<view_type_a_ll, Device>(1024);
  // Test::impl_test_nrm2<view_type_a_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA*, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrm2w<view_type_a_lr, Device>(0);
  Test::impl_test_nrm2w<view_type_a_lr, Device>(13);
  Test::impl_test_nrm2w<view_type_a_lr, Device>(1024);
  // Test::impl_test_nrm2<view_type_a_lr, Device>(132231);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
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
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_nrm2w_mv<view_type_a_ll, Device>(0, 5);
  Test::impl_test_nrm2w_mv<view_type_a_ll, Device>(13, 5);
  Test::impl_test_nrm2w_mv<view_type_a_ll, Device>(1024, 5);
  Test::impl_test_nrm2w_mv<view_type_a_ll, Device>(789, 1);
  // Test::impl_test_nrm2w_mv<view_type_a_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_nrm2w_mv<view_type_a_lr, Device>(0, 5);
  Test::impl_test_nrm2w_mv<view_type_a_lr, Device>(13, 5);
  Test::impl_test_nrm2w_mv<view_type_a_lr, Device>(1024, 5);
  Test::impl_test_nrm2w_mv<view_type_a_lr, Device>(789, 1);
  // Test::impl_test_nrm2w_mv<view_type_a_lr, Device>(132231,5);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
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
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm2w_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_float");
  test_nrm2w<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm2w_mv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_mv_float");
  test_nrm2w_mv<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm2w_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_double");
  test_nrm2w<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm2w_mv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_mv_double");
  test_nrm2w_mv<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm2w_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_complex_double");
  test_nrm2w<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm2w_mv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_mv_complex_double");
  test_nrm2w_mv<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, nrm2w_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_int");
  test_nrm2w<int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
TEST_F(TestCategory, nrm2w_mv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::nrm2w_mv_int");
  test_nrm2w_mv<int, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
