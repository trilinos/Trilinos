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
#include <KokkosBlas1_team_nrm2.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ViewTypeA, class Device>
void impl_test_team_nrm2(int N, int K) {
  using execution_space = typename Device::execution_space;
  typedef Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;

  // Launch K teams of the maximum number of threads per team
  const team_policy policy(K, Kokkos::AUTO);

  typedef typename ViewTypeA::value_type ScalarA;
  typedef Kokkos::ArithTraits<ScalarA> AT;

  view_stride_adapter<ViewTypeA> a("A", N, K);

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

  Kokkos::fill_random(a.d_view, rand_pool, ScalarA(10));

  Kokkos::deep_copy(a.h_base, a.d_base);

  typename AT::mag_type *expected_result = new typename AT::mag_type[K];
  for (int j = 0; j < K; j++) {
    expected_result[j] = typename AT::mag_type();
    for (int i = 0; i < N; i++) expected_result[j] += AT::abs(a.h_view(i, j)) * AT::abs(a.h_view(i, j));
    expected_result[j] = Kokkos::ArithTraits<typename AT::mag_type>::sqrt(expected_result[j]);
  }

  double eps = std::is_same<ScalarA, float>::value ? 2 * 1e-5 : 1e-7;

  Kokkos::View<typename AT::mag_type *, Kokkos::HostSpace> r("Nrm2::Result", K);
  Kokkos::View<typename AT::mag_type *, Device> d_r("Nrm2::Result", K);

  // KokkosBlas::nrm2(r,a);
  Kokkos::parallel_for(
      "KokkosBlas::Test::TeamNrm2", policy, KOKKOS_LAMBDA(const team_member &teamMember) {
        const int teamId = teamMember.league_rank();
        d_r(teamId)      = KokkosBlas::Experimental::nrm2(teamMember, Kokkos::subview(a.d_view, Kokkos::ALL(), teamId));
      });
  Kokkos::deep_copy(r, d_r);
  for (int k = 0; k < K; k++) {
    typename AT::mag_type nonconst_result = r(k);
    EXPECT_NEAR_KK(nonconst_result, expected_result[k], eps * expected_result[k]);
  }

  // KokkosBlas::nrm2(r,c_a);
  Kokkos::parallel_for(
      "KokkosBlas::Test::TeamNrm2", policy, KOKKOS_LAMBDA(const team_member &teamMember) {
        const int teamId = teamMember.league_rank();
        d_r(teamId) =
            KokkosBlas::Experimental::nrm2(teamMember, Kokkos::subview(a.d_view_const, Kokkos::ALL(), teamId));
      });
  Kokkos::deep_copy(r, d_r);
  for (int k = 0; k < K; k++) {
    typename AT::mag_type const_result = r(k);
    EXPECT_NEAR_KK(const_result, expected_result[k], eps * expected_result[k]);
  }

  delete[] expected_result;
}
}  // namespace Test

template <class ScalarA, class Device>
int test_team_nrm2() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA **, Kokkos::LayoutLeft, Device> view_type_a_ll;
  Test::impl_test_team_nrm2<view_type_a_ll, Device>(0, 5);
  Test::impl_test_team_nrm2<view_type_a_ll, Device>(13, 5);
  Test::impl_test_team_nrm2<view_type_a_ll, Device>(124, 5);
  // Test::impl_test_team_nrm2<view_type_a_ll, Device>(132231,5);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA **, Kokkos::LayoutRight, Device> view_type_a_lr;
  Test::impl_test_team_nrm2<view_type_a_lr, Device>(0, 5);
  Test::impl_test_team_nrm2<view_type_a_lr, Device>(13, 5);
  Test::impl_test_team_nrm2<view_type_a_lr, Device>(124, 5);
  // Test::impl_test_team_nrm2<view_type_a_lr, Device>(132231,5);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA **, Kokkos::LayoutStride, Device> view_type_a_ls;
  Test::impl_test_team_nrm2<view_type_a_ls, Device>(0, 5);
  Test::impl_test_team_nrm2<view_type_a_ls, Device>(13, 5);
  Test::impl_test_team_nrm2<view_type_a_ls, Device>(124, 5);
  // Test::impl_test_team_nrm2<view_type_a_ls, Device>(132231,5);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_nrm2_float) { test_team_nrm2<float, TestDevice>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_nrm2_double) { test_team_nrm2<double, TestDevice>(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_nrm2_complex_double) { test_team_nrm2<Kokkos::complex<double>, TestDevice>(); }
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, team_nrm2_int) { test_team_nrm2<int, TestDevice>(); }
#endif

#endif  // Check for lambda availability in CUDA backend
