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
#include <KokkosBlas1_team_mult.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class Device>
void impl_test_team_mult(int N) {
  using execution_space = typename Device::execution_space;
  typedef Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;

  // Launch M teams of the maximum number of threads per team
  int M = 4;
  const team_policy policy(M, Kokkos::AUTO);
  const int team_data_siz = (N % M == 0) ? (N / M) : (N / M + 1);

  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef typename ViewTypeC::value_type ScalarC;

  ScalarA a  = 3;
  ScalarB b  = 5;
  double eps = std::is_same<ScalarC, float>::value ? 2 * 1e-5 : 1e-7;

  view_stride_adapter<ViewTypeA> x("X", N);
  view_stride_adapter<ViewTypeB> y("Y", N);
  view_stride_adapter<ViewTypeC> z("Z", N);
  view_stride_adapter<ViewTypeC> org_z("Org_Z", N);

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

  Kokkos::fill_random(x.d_view, rand_pool, ScalarA(10));
  Kokkos::fill_random(y.d_view, rand_pool, ScalarB(10));
  Kokkos::fill_random(z.d_view, rand_pool, ScalarC(10));

  Kokkos::deep_copy(org_z.h_base, z.d_base);

  Kokkos::deep_copy(x.h_base, x.d_base);
  Kokkos::deep_copy(y.h_base, y.d_base);
  Kokkos::deep_copy(z.h_base, z.d_base);

  ScalarA expected_result = 0;
  for (int i = 0; i < N; i++)
    expected_result += ScalarC(b * z.h_view(i) + a * x.h_view(i) * y.h_view(i)) *
                       ScalarC(b * z.h_view(i) + a * x.h_view(i) * y.h_view(i));

  // KokkosBlas::mult(b,z,a,x,y);
  Kokkos::parallel_for(
      "KokkosBlas::Test::TeamMult", policy, KOKKOS_LAMBDA(const team_member &teamMember) {
        const int teamId = teamMember.league_rank();
        KokkosBlas::Experimental::mult(
            teamMember, b,
            Kokkos::subview(z.d_view, Kokkos::make_pair(teamId * team_data_siz,
                                                        (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)),
            a,
            Kokkos::subview(x.d_view, Kokkos::make_pair(teamId * team_data_siz,
                                                        (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)),
            Kokkos::subview(y.d_view, Kokkos::make_pair(teamId * team_data_siz,
                                                        (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)));
      });
  ScalarC nonconst_nonconst_result = KokkosBlas::dot(z.d_view, z.d_view);
  EXPECT_NEAR_KK(nonconst_nonconst_result, expected_result, eps * expected_result);

  // Reset z on device to orig and run again with const-valued y
  Kokkos::deep_copy(z.d_base, org_z.h_base);
  // KokkosBlas::mult(b,z,a,x,c_y);
  Kokkos::parallel_for(
      "KokkosBlas::Test::TeamMult", policy, KOKKOS_LAMBDA(const team_member &teamMember) {
        const int teamId = teamMember.league_rank();
        KokkosBlas::Experimental::mult(
            teamMember, b,
            Kokkos::subview(z.d_view, Kokkos::make_pair(teamId * team_data_siz,
                                                        (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)),
            a,
            Kokkos::subview(x.d_view, Kokkos::make_pair(teamId * team_data_siz,
                                                        (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)),
            Kokkos::subview(y.d_view_const, Kokkos::make_pair(teamId * team_data_siz,
                                                              (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)));
      });
  ScalarC const_nonconst_result = KokkosBlas::dot(z.d_view, z.d_view);
  EXPECT_NEAR_KK(const_nonconst_result, expected_result, eps * expected_result);

  // Reset z again to orig, and run with both x and y const
  Kokkos::deep_copy(z.d_base, org_z.h_base);
  // KokkosBlas::mult(b,z,a,c_x,c_y);
  Kokkos::parallel_for(
      "KokkosBlas::Test::TeamMult", policy, KOKKOS_LAMBDA(const team_member &teamMember) {
        const int teamId = teamMember.league_rank();
        KokkosBlas::Experimental::mult(
            teamMember, b,
            Kokkos::subview(z.d_view, Kokkos::make_pair(teamId * team_data_siz,
                                                        (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)),
            a,
            Kokkos::subview(x.d_view_const, Kokkos::make_pair(teamId * team_data_siz,
                                                              (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)),
            Kokkos::subview(y.d_view_const, Kokkos::make_pair(teamId * team_data_siz,
                                                              (teamId < M - 1) ? (teamId + 1) * team_data_siz : N)));
      });
  ScalarC const_const_result = KokkosBlas::dot(z.d_view, z.d_view);
  EXPECT_NEAR_KK(const_const_result, expected_result, eps * expected_result);
}

template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class Device>
void impl_test_team_mult_mv(int N, int K) {
  using execution_space = typename Device::execution_space;
  typedef Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;

  // Launch K teams of the maximum number of threads per team
  const team_policy policy(K, Kokkos::AUTO);

  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef typename ViewTypeC::value_type ScalarC;

  // x is rank-1, all others are rank-2
  view_stride_adapter<ViewTypeA> x("X", N);
  view_stride_adapter<ViewTypeB> y("Y", N, K);
  view_stride_adapter<ViewTypeC> z("Z", N, K);
  view_stride_adapter<ViewTypeC> org_z("Org_Z", N, K);

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

  typename Kokkos::ArithTraits<ScalarC>::mag_type const max_val = 10;
  Kokkos::fill_random(x.d_view, rand_pool, ScalarA(max_val));
  Kokkos::fill_random(y.d_view, rand_pool, ScalarB(max_val));
  Kokkos::fill_random(z.d_view, rand_pool, ScalarC(max_val));

  Kokkos::deep_copy(org_z.h_base, z.d_base);

  Kokkos::deep_copy(x.h_base, x.d_base);
  Kokkos::deep_copy(y.h_base, y.d_base);

  ScalarA a = 3;
  ScalarB b = 5;

  // In the operation z = (b*z) + (a*x*y) we estimate
  // the largest rounding error to be dominated by max(b*z, a*x*y)
  // Since b and a are known and the largest value in z, x and y
  // is set by the variables max_val, the error upper bound will be
  //         max_error = a * max_val * max_val
  typename Kokkos::ArithTraits<ScalarC>::mag_type const eps = Kokkos::ArithTraits<ScalarC>::epsilon();
  typename Kokkos::ArithTraits<ScalarC>::mag_type const max_error =
      Kokkos::ArithTraits<ScalarC>::abs(a) * max_val * max_val * eps;

  // KokkosBlas::mult(b,z,a,x,y);
  Kokkos::parallel_for(
      "KokkosBlas::Test::TeamMult", policy, KOKKOS_LAMBDA(const team_member &teamMember) {
        const int teamId = teamMember.league_rank();
        KokkosBlas::Experimental::mult(teamMember, b, Kokkos::subview(z.d_view, Kokkos::ALL(), teamId), a, x.d_view,
                                       Kokkos::subview(y.d_view, Kokkos::ALL(), teamId));
      });

  Kokkos::deep_copy(z.h_base, z.d_base);

  ScalarC temp;
  for (int j = 0; j < K; j++) {
    for (int i = 0; i < N; i++) {
      temp = ScalarC(b * org_z.h_view(i, j) + a * x.h_view(i) * y.h_view(i, j));
      EXPECT_NEAR_KK(temp, z.h_view(i, j), max_error);
    }
  }

  // Reset z on device and run again with const y
  Kokkos::deep_copy(z.d_base, org_z.h_base);
  // KokkosBlas::mult(b,z,a,x,c_y);
  Kokkos::parallel_for(
      "KokkosBlas::Test::TeamMult", policy, KOKKOS_LAMBDA(const team_member &teamMember) {
        const int teamId = teamMember.league_rank();
        KokkosBlas::Experimental::mult(teamMember, b, Kokkos::subview(z.d_view, Kokkos::ALL(), teamId), a, x.d_view,
                                       Kokkos::subview(y.d_view_const, Kokkos::ALL(), teamId));
      });
  Kokkos::deep_copy(z.h_base, z.d_base);

  for (int k = 0; k < K; k++) {
    for (int i = 0; i < N; ++i) {
      temp = ScalarC(b * org_z.h_view(i, k) + a * x.h_view(i) * y.h_view(i, k));
      EXPECT_NEAR_KK(temp, z.h_view(i, k), max_error);
    }
  }
}
}  // namespace Test

template <class ScalarA, class ScalarB, class ScalarC, class Device>
int test_team_mult() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA *, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB *, Kokkos::LayoutLeft, Device> view_type_b_ll;
  typedef Kokkos::View<ScalarC *, Kokkos::LayoutLeft, Device> view_type_c_ll;
  Test::impl_test_team_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(0);
  Test::impl_test_team_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(13);
  Test::impl_test_team_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(124);
  // Test::impl_test_team_mult<view_type_a_ll, view_type_b_ll, view_type_c_ll,
  // Device>(132231);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA *, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB *, Kokkos::LayoutRight, Device> view_type_b_lr;
  typedef Kokkos::View<ScalarC *, Kokkos::LayoutRight, Device> view_type_c_lr;
  Test::impl_test_team_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(0);
  Test::impl_test_team_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(13);
  Test::impl_test_team_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(124);
  // Test::impl_test_team_mult<view_type_a_lr, view_type_b_lr, view_type_c_lr,
  // Device>(132231);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA *, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB *, Kokkos::LayoutStride, Device> view_type_b_ls;
  typedef Kokkos::View<ScalarC *, Kokkos::LayoutStride, Device> view_type_c_ls;
  Test::impl_test_team_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(0);
  Test::impl_test_team_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(13);
  Test::impl_test_team_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(124);
  // Test::impl_test_team_mult<view_type_a_ls, view_type_b_ls, view_type_c_ls,
  // Device>(132231);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_team_mult<view_type_a_ls, view_type_b_ll, view_type_c_lr, Device>(124);
  Test::impl_test_team_mult<view_type_a_ll, view_type_b_ls, view_type_c_lr, Device>(124);
#endif

  return 1;
}

template <class ScalarA, class ScalarB, class ScalarC, class Device>
int test_team_mult_mv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA *, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarB **, Kokkos::LayoutLeft, Device> view_type_b_ll;
  typedef Kokkos::View<ScalarC **, Kokkos::LayoutLeft, Device> view_type_c_ll;
  Test::impl_test_team_mult_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(0, 5);
  Test::impl_test_team_mult_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(13, 5);
  Test::impl_test_team_mult_mv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(124, 5);
  // Test::impl_test_team_mult_mv<view_type_a_ll, view_type_b_ll,
  // view_type_c_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA *, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarB **, Kokkos::LayoutRight, Device> view_type_b_lr;
  typedef Kokkos::View<ScalarC **, Kokkos::LayoutRight, Device> view_type_c_lr;
  Test::impl_test_team_mult_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(0, 5);
  Test::impl_test_team_mult_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(13, 5);
  Test::impl_test_team_mult_mv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(124, 5);
  // Test::impl_test_team_mult_mv<view_type_a_lr, view_type_b_lr,
  // view_type_c_lr, Device>(132231,5);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA *, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarB **, Kokkos::LayoutStride, Device> view_type_b_ls;
  typedef Kokkos::View<ScalarC **, Kokkos::LayoutStride, Device> view_type_c_ls;
  Test::impl_test_team_mult_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(0, 5);
  Test::impl_test_team_mult_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(13, 5);
  Test::impl_test_team_mult_mv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(124, 5);
  // Test::impl_test_team_mult_mv<view_type_a_ls, view_type_b_ls,
  // view_type_c_ls, Device>(132231,5);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_team_mult_mv<view_type_a_ls, view_type_b_ll, view_type_c_lr, Device>(124, 5);
  Test::impl_test_team_mult_mv<view_type_a_ll, view_type_b_ls, view_type_c_lr, Device>(124, 5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_mult_float) { test_team_mult<float, float, float, TestDevice>(); }
TEST_F(TestCategory, team_mult_mv_float) { test_team_mult_mv<float, float, float, TestDevice>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_mult_double) { test_team_mult<double, double, double, TestDevice>(); }
TEST_F(TestCategory, team_mult_mv_double) { test_team_mult_mv<double, double, double, TestDevice>(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_mult_complex_double) {
  test_team_mult<Kokkos::complex<double>, Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
}
TEST_F(TestCategory, team_mult_mv_complex_double) {
  test_team_mult_mv<Kokkos::complex<double>, Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_mult_int) { test_team_mult<int, int, int, TestDevice>(); }
TEST_F(TestCategory, team_mult_mv_int) { test_team_mult_mv<int, int, int, TestDevice>(); }
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, team_mult_double_int) { test_team_mult<double, int, float, TestDevice>(); }
TEST_F(TestCategory, team_mult_double_mv_int) { test_team_mult_mv<double, int, float, TestDevice>(); }
#endif

#endif  // Check for lambda availability in CUDA backend
