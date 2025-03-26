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

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_trsv_l_nt_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_nt_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_t_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_t_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_c_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_c_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}

TEST_F(TestCategory, batched_scalar_serial_trsv_u_nt_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_nt_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_t_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_t_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_c_u_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_c_n_fcomplex_fcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_trsv_l_nt_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_nt_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_t_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_t_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_c_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_l_c_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Lower, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}

TEST_F(TestCategory, batched_scalar_serial_trsv_u_nt_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_nt_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::NoTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_t_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_t_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::Transpose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_c_u_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::ConjTranspose, Diag::Unit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_trsv_u_c_n_dcomplex_dcomplex) {
  using param_tag_type = ::Test::Trsv::ParamTag<Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit>;
  using algo_tag_type  = Algo::Trsv::Unblocked;
  test_batched_trsv<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
#endif
