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

// We do not ETI half-types. Only test this if ETI ONLY is off
// and bhalf_t is not an alias to float.
#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS) && \
    defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
/********************* BatchLayout::Left *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_bhalf_bhalf_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_bhalf_bhalf_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_bhalf_bhalf_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_bhalf_bhalf_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type>();
}
/********************* BatchLayout::Right *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_bhalf_bhalf_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_bhalf_bhalf_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_bhalf_bhalf_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_bhalf_bhalf_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type>();
}
#endif  // KOKKOS_BHALF_T_IS_FLOAT

// We do not ETI half-types. Only test this if ETI ONLY is off
// and half_t is not an alias to float.
#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS) && \
    defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
/********************* BatchLayout::Left *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_half_half_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_half_half_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_half_half_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_half_half_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type>();
}
/********************* BatchLayout::Right *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_half_half_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_half_half_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_half_half_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_half_half_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type>();
}
#endif  // KOKKOS_HALF_T_IS_FLOAT

#if defined(KOKKOSKERNELS_INST_FLOAT)
/********************* BatchLayout::Left *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_float_float_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_float_float_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_float_float_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_float_float_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, float, float, param_tag_type>();
}
/********************* BatchLayout::Right *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_float_float_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_float_float_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_float_float_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, float, float, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_float_float_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, float, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
/********************* BatchLayout::Left *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_double_double_left) {
  using param_tag_type = ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Left>;

  test_batched_gemm<TestDevice, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_double_double_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_double_double_left) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_double_double_left) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Left> param_tag_type;

  test_batched_gemm<TestDevice, double, double, param_tag_type>();
}
/********************* BatchLayout::Right *********************/
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_nt_double_double_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_nt_double_double_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::NoTranspose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_nt_t_double_double_right) {
  typedef ::Test::SharedParamTag<Trans::NoTranspose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, double, double, param_tag_type>();
}
TEST_F(TestCategory, batched_scalar_batched_gemm_t_t_double_double_right) {
  typedef ::Test::SharedParamTag<Trans::Transpose, Trans::Transpose, BatchLayout::Right> param_tag_type;

  test_batched_gemm<TestDevice, double, double, param_tag_type>();
}
#endif
