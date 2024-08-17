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
// Note: Luc Berger-Vergiat 04/14/21
//       This tests uses KOKKOS_LAMBDA so we need
//       to make sure that these are enabled in
//       the CUDA backend before including this test.
#if !defined(TEST_CUDA_BLAS_CPP) || defined(KOKKOS_ENABLE_CUDA_LAMBDA)

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas1_team_axpby.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_team_axpy(int N) {
  using execution_space = typename Device::execution_space;
  typedef Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;

  // Launch M teams of the maximum number of threads per team
  int M = 4;
  const team_policy policy(M, Kokkos::AUTO);
  const int team_data_siz = (N % M == 0) ? (N / M) : (N / M + 1);

  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;

  view_stride_adapter<ViewTypeA> x("X", N);
  view_stride_adapter<ViewTypeB> y("Y", N);
  view_stride_adapter<ViewTypeB> org_y("Y", N);

  ScalarA a  = 3;
  double eps = std::is_same<ScalarA, float>::value ? 2 * 1e-5 : 1e-7;

  Kokkos::Random_XorShift64_Pool<typename Device::execution_space> rand_pool(13718);

  Kokkos::fill_random(x.d_view, rand_pool, ScalarA(10));
  Kokkos::fill_random(y.d_view, rand_pool, ScalarB(10));

  Kokkos::deep_copy(x.h_base, x.d_base);
  Kokkos::deep_copy(y.h_base, y.d_base);
  Kokkos::deep_copy(org_y.h_base, y.d_base);

  ScalarA expected_result = 0;
  for (int i = 0; i < N; i++)
    expected_result += ScalarB(a * x.h_view(i) + y.h_view(i)) * ScalarB(a * x.h_view(i) + y.h_view(i));

  // KokkosBlas::axpy(a,x,y);
  Kokkos::parallel_for(
      "KokkosBlas::Test::TeamAxpy", policy, KOKKOS_LAMBDA(const team_member &teamMember) {
        const int teamId = teamMember.league_rank();
        KokkosBlas::Experimental::axpy(
            teamMember, a,
            Kokkos::subview(x.d_view, Kokkos::make_pair(teamId * team_data_siz,
                                                        (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)),
            Kokkos::subview(y.d_view, Kokkos::make_pair(teamId * team_data_siz,
                                                        (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)));
      });

  ScalarB nonconst_nonconst_result = KokkosBlas::dot(y.d_view, y.d_view);
  EXPECT_NEAR_KK(nonconst_nonconst_result, expected_result, eps * expected_result);

  Kokkos::deep_copy(y.d_base, org_y.h_base);

  // KokkosBlas::axpy(a,c_x,y);
  Kokkos::parallel_for(
      "KokkosBlas::Test::TeamAxpy", policy, KOKKOS_LAMBDA(const team_member &teamMember) {
        const int teamId = teamMember.league_rank();
        KokkosBlas::Experimental::axpy(
            teamMember, a,
            Kokkos::subview(x.d_view_const, Kokkos::make_pair(teamId * team_data_siz,
                                                              (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)),
            Kokkos::subview(y.d_view, Kokkos::make_pair(teamId * team_data_siz,
                                                        (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)));
      });

  ScalarB const_nonconst_result = KokkosBlas::dot(y.d_view_const, y.d_view_const);
  EXPECT_NEAR_KK(const_nonconst_result, expected_result, eps * expected_result);
}

template <class ViewTypeA, class ViewTypeB, class Device>
void impl_test_team_axpy_mv(int N, int K) {
  using execution_space = typename Device::execution_space;
  typedef Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;

  // Launch K teams of the maximum number of threads per team
  const team_policy policy(K, Kokkos::AUTO);

  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;

  view_stride_adapter<ViewTypeA> x("X", N, K);
  view_stride_adapter<ViewTypeB> y("Y", N, K);
  view_stride_adapter<ViewTypeB> org_y("Org_Y", N, K);

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

  Kokkos::fill_random(x.d_view, rand_pool, ScalarA(10));
  Kokkos::fill_random(y.d_view, rand_pool, ScalarB(10));

  Kokkos::deep_copy(x.h_base, x.d_base);
  Kokkos::deep_copy(y.h_base, y.d_base);
  Kokkos::deep_copy(org_y.h_base, y.d_base);

  ScalarA a = 3;

  ScalarA *expected_result = new ScalarA[K];
  for (int j = 0; j < K; j++) {
    expected_result[j] = ScalarA();
    for (int i = 0; i < N; i++)
      expected_result[j] += ScalarB(a * x.h_view(i, j) + y.h_view(i, j)) * ScalarB(a * x.h_view(i, j) + y.h_view(i, j));
  }

  double eps = std::is_same<ScalarA, float>::value ? 2 * 1e-5 : 1e-7;

  Kokkos::View<ScalarB *, Kokkos::HostSpace> r("Dot::Result", K);

  // KokkosBlas::axpy(a,x,y);
  Kokkos::parallel_for(
      "KokkosBlas::Test::TeamAxpy", policy, KOKKOS_LAMBDA(const team_member &teamMember) {
        const int teamId = teamMember.league_rank();
        KokkosBlas::Experimental::axpy(teamMember, a, Kokkos::subview(x.d_view, Kokkos::ALL(), teamId),
                                       Kokkos::subview(y.d_view, Kokkos::ALL(), teamId));
      });

  KokkosBlas::dot(r, y.d_view, y.d_view);
  for (int k = 0; k < K; k++) {
    ScalarA nonconst_nonconst_result = r(k);
    EXPECT_NEAR_KK(nonconst_nonconst_result, expected_result[k], eps * expected_result[k]);
  }

  Kokkos::deep_copy(y.d_base, org_y.h_base);

  // KokkosBlas::axpy(a,c_x,y);
  Kokkos::parallel_for(
      "KokkosBlas::Test::TeamAxpy", policy, KOKKOS_LAMBDA(const team_member &teamMember) {
        const int teamId = teamMember.league_rank();
        KokkosBlas::Experimental::axpy(teamMember, a, Kokkos::subview(x.d_view_const, Kokkos::ALL(), teamId),
                                       Kokkos::subview(y.d_view, Kokkos::ALL(), teamId));
      });

  KokkosBlas::dot(r, y.d_view, y.d_view);
  for (int k = 0; k < K; k++) {
    ScalarA const_non_const_result = r(k);
    EXPECT_NEAR_KK(const_non_const_result, expected_result[k], eps * expected_result[k]);
  }

  delete[] expected_result;
}
}  // namespace Test

template <class ScalarA, class ScalarB, class Device>
int test_team_axpy() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA *, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB *, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_team_axpy<view_type_a_ll, view_type_b_ll, Device>(0);
  Test::impl_test_team_axpy<view_type_a_ll, view_type_b_ll, Device>(13);
  Test::impl_test_team_axpy<view_type_a_ll, view_type_b_ll, Device>(124);
  // Test::impl_test_team_axpy<view_type_a_ll, view_type_b_ll, Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA *, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB *, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_team_axpy<view_type_a_lr, view_type_b_lr, Device>(0);
  Test::impl_test_team_axpy<view_type_a_lr, view_type_b_lr, Device>(13);
  Test::impl_test_team_axpy<view_type_a_lr, view_type_b_lr, Device>(124);
  // Test::impl_test_team_axpy<view_type_a_lr, view_type_b_lr, Device>(132231);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA *, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB *, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_team_axpy<view_type_a_ls, view_type_b_ls, Device>(0);
  Test::impl_test_team_axpy<view_type_a_ls, view_type_b_ls, Device>(13);
  Test::impl_test_team_axpy<view_type_a_ls, view_type_b_ls, Device>(124);
  // Test::impl_test_team_axpy<view_type_a_ls, view_type_b_ls, Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_team_axpy<view_type_a_ls, view_type_b_ll, Device>(124);
  Test::impl_test_team_axpy<view_type_a_ll, view_type_b_ls, Device>(124);
#endif

  return 1;
}

template <class ScalarA, class ScalarB, class Device>
int test_team_axpy_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA **, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB **, Kokkos::LayoutLeft, Device> view_type_b_ll;
  Test::impl_test_team_axpy_mv<view_type_a_ll, view_type_b_ll, Device>(0, 5);
  Test::impl_test_team_axpy_mv<view_type_a_ll, view_type_b_ll, Device>(13, 5);
  Test::impl_test_team_axpy_mv<view_type_a_ll, view_type_b_ll, Device>(124, 5);
  // Test::impl_test_team_axpy_mv<view_type_a_ll, view_type_b_ll,
  // Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA **, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB **, Kokkos::LayoutRight, Device> view_type_b_lr;
  Test::impl_test_team_axpy_mv<view_type_a_lr, view_type_b_lr, Device>(0, 5);
  Test::impl_test_team_axpy_mv<view_type_a_lr, view_type_b_lr, Device>(13, 5);
  Test::impl_test_team_axpy_mv<view_type_a_lr, view_type_b_lr, Device>(124, 5);
  // Test::impl_test_team_axpy_mv<view_type_a_lr, view_type_b_lr,
  // Device>(132231,5);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA **, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB **, Kokkos::LayoutStride, Device> view_type_b_ls;
  Test::impl_test_team_axpy_mv<view_type_a_ls, view_type_b_ls, Device>(0, 5);
  Test::impl_test_team_axpy_mv<view_type_a_ls, view_type_b_ls, Device>(13, 5);
  Test::impl_test_team_axpy_mv<view_type_a_ls, view_type_b_ls, Device>(124, 5);
  // Test::impl_test_team_axpy_mv<view_type_a_ls, view_type_b_ls,
  // Device>(132231,5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_team_axpy_mv<view_type_a_ls, view_type_b_ll, Device>(124, 5);
  Test::impl_test_team_axpy_mv<view_type_a_ll, view_type_b_ls, Device>(124, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_axpy_float) { test_team_axpy<float, float, TestDevice>(); }
TEST_F(TestCategory, team_axpy_mv_float) { test_team_axpy_mv<float, float, TestDevice>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_axpy_double) { test_team_axpy<double, double, TestDevice>(); }
TEST_F(TestCategory, team_axpy_mv_double) { test_team_axpy_mv<double, double, TestDevice>(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_axpy_complex_double) {
  test_team_axpy<Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
}
TEST_F(TestCategory, team_axpy_mv_complex_double) {
  test_team_axpy_mv<Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_axpy_int) { test_team_axpy<int, int, TestDevice>(); }
TEST_F(TestCategory, team_axpy_mv_int) { test_team_axpy_mv<int, int, TestDevice>(); }
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, team_axpy_double_int) { test_team_axpy<double, int, TestDevice>(); }
TEST_F(TestCategory, team_axpy_double_mv_int) { test_team_axpy_mv<double, int, TestDevice>(); }
#endif

#endif  // Check for lambda availability in CUDA backend
