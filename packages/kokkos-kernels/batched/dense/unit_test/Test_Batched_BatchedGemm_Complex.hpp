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
/********************* BatchLayout::Left *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_scomplex_scomplex_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_scomplex_scomplex_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_scomplex_scomplex_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_scomplex_scomplex_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type>();
}
/********************* BatchLayout::Right *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_scomplex_scomplex_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_scomplex_scomplex_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_scomplex_scomplex_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_scomplex_scomplex_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
/********************* BatchLayout::Left *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_dcomplex_dcomplex_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_dcomplex_dcomplex_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_dcomplex_dcomplex_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_dcomplex_dcomplex_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type>();
}
/********************* BatchLayout::Right *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_dcomplex_dcomplex_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_dcomplex_dcomplex_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_dcomplex_dcomplex_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_dcomplex_dcomplex_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type>();
}
#endif
