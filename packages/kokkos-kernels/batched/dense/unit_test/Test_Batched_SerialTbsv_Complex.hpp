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
// NO TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_nt_u_fcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_nt_n_fcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_u_fcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_n_fcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
// TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_t_u_fcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_t_n_fcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::NonUnit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_u_fcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_n_fcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::NonUnit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}

// CONJUGATE TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_c_u_fcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_c_n_fcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_c_u_fcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_c_n_fcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// NO TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_nt_u_dcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_nt_n_dcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_u_dcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_nt_n_dcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
// TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_t_u_dcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_t_n_dcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::NonUnit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_u_dcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::Unit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_t_n_dcomplex) {
  using param_tag_type =
      ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::NonUnit>;
  using algo_tag_type = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}

// CONJUGATE TRANSPOSE
TEST_F(TestCategory, batched_serial_tbsv_l_c_u_dcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_l_c_n_dcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_c_u_dcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Diag::Unit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_serial_tbsv_u_c_n_dcomplex) {
  using param_tag_type = ::Test::Tbsv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Diag::NonUnit>;
  using algo_tag_type  = typename Algo::Tbsv::Unblocked;

  test_batched_tbsv<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
#endif
