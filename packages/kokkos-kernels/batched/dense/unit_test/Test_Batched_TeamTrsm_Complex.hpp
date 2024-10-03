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

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, batched_scalar_team_trsm_l_l_nt_u_dcomplex_dcomplex) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_l_nt_n_dcomplex_dcomplex) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_u_nt_u_dcomplex_dcomplex) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_u_nt_n_dcomplex_dcomplex) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_r_u_nt_u_dcomplex_dcomplex) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_r_u_nt_n_dcomplex_dcomplex) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
//
TEST_F(TestCategory, batched_scalar_team_trsm_l_l_t_u_dcomplex_dcomplex) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_l_t_n_dcomplex_dcomplex) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_u_t_u_dcomplex_dcomplex) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_u_t_n_dcomplex_dcomplex) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}

TEST_F(TestCategory, batched_scalar_team_trsm_l_l_nt_u_dcomplex_double) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_l_nt_n_dcomplex_double) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_u_nt_u_dcomplex_double) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_u_nt_n_dcomplex_double) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_r_u_nt_u_dcomplex_double) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_r_u_nt_n_dcomplex_double) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
//
TEST_F(TestCategory, batched_scalar_team_trsm_l_l_t_u_dcomplex_double) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_l_t_n_dcomplex_double) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Lower, Trans::Transpose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_u_t_u_dcomplex_double) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::Unit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_team_trsm_l_u_t_n_dcomplex_double) {
  typedef ::Test::TeamTrsm::ParamTag<Side::Left, Uplo::Upper, Trans::Transpose, Diag::NonUnit> param_tag_type;
  typedef Algo::Trsm::Blocked algo_tag_type;
  test_batched_team_trsm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
#endif
