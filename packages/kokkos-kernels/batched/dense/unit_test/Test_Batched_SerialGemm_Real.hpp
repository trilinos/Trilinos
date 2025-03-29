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
#if defined(KOKKOS_BHALF_T_IS_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_nt_bhalf_bhalf) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_nt_bhalf_bhalf) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_t_bhalf_bhalf) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_t_bhalf_bhalf) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;

  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, ::Test::bhalfScalarType, ::Test::bhalfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
#endif  // KOKKOS_BHALF_T_IS_FLOAT

#if defined(KOKKOS_HALF_T_IS_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_nt_half_half) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_nt_half_half) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_t_half_half) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_t_half_half) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;

  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, ::Test::halfScalarType, ::Test::halfScalarType, param_tag_type,
                    KokkosBatched::Algo::Gemm::Unblocked>();
}
#endif  // KOKKOS_HALF_T_IS_FLOAT

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_nt_float_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_nt_float_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_t_float_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_t_float_float) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, float, float, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_nt_double_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}

TEST_F(TestCategory, batched_scalar_serial_gemm_t_nt_double_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose>;
  test_batched_gemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_t_double_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_t_double_double) {
  using param_tag_type = ::Test::Gemm::ParamTag<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::Transpose>;
  test_batched_gemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Blocked>();
  test_batched_gemm<TestDevice, double, double, param_tag_type, KokkosBatched::Algo::Gemm::Unblocked>();
}
#endif
