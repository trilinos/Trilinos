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

/// dcomplex, dcomplex

TEST_F(TestCategory, batched_scalar_serial_gemm_nt_nt_dcomplex_dcomplex) {
  typedef ::Test::Gemm::ParamTag<Trans::NoTranspose, Trans::NoTranspose> param_tag_type;
  typedef Algo::Gemm::Blocked algo_tag_type;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_nt_dcomplex_dcomplex) {
  typedef ::Test::Gemm::ParamTag<Trans::Transpose, Trans::NoTranspose> param_tag_type;
  typedef Algo::Gemm::Blocked algo_tag_type;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_t_dcomplex_dcomplex) {
  typedef ::Test::Gemm::ParamTag<Trans::NoTranspose, Trans::Transpose> param_tag_type;
  typedef Algo::Gemm::Blocked algo_tag_type;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_t_dcomplex_dcomplex) {
  typedef ::Test::Gemm::ParamTag<Trans::Transpose, Trans::Transpose> param_tag_type;
  typedef Algo::Gemm::Blocked algo_tag_type;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
// TEST_F( TestCategory, batched_scalar_serial_gemm_ct_nt_dcomplex_dcomplex ) {
//   typedef ::Test::Gemm::ParamTag<Trans::ConjTranspose,Trans::NoTranspose>
//   param_tag_type; typedef Algo::Gemm::Blocked algo_tag_type;
//   test_batched_gemm<TestDevice,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>();
// }
// TEST_F( TestCategory, batched_scalar_serial_gemm_nt_ct_dcomplex_dcomplex ) {
//   typedef ::Test::Gemm::ParamTag<Trans::NoTranspose,Trans::ConjTranspose>
//   param_tag_type; typedef Algo::Gemm::Blocked algo_tag_type;
//   test_batched_gemm<TestDevice,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,algo_tag_type>();
// }

/// dcomplex, double

TEST_F(TestCategory, batched_scalar_serial_gemm_nt_nt_dcomplex_double) {
  typedef ::Test::Gemm::ParamTag<Trans::NoTranspose, Trans::NoTranspose> param_tag_type;
  typedef Algo::Gemm::Blocked algo_tag_type;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_nt_dcomplex_double) {
  typedef ::Test::Gemm::ParamTag<Trans::Transpose, Trans::NoTranspose> param_tag_type;
  typedef Algo::Gemm::Blocked algo_tag_type;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_nt_t_dcomplex_double) {
  typedef ::Test::Gemm::ParamTag<Trans::NoTranspose, Trans::Transpose> param_tag_type;
  typedef Algo::Gemm::Blocked algo_tag_type;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, batched_scalar_serial_gemm_t_t_dcomplex_double) {
  typedef ::Test::Gemm::ParamTag<Trans::Transpose, Trans::Transpose> param_tag_type;
  typedef Algo::Gemm::Blocked algo_tag_type;
  test_batched_gemm<TestDevice, Kokkos::complex<double>, double, param_tag_type, algo_tag_type>();
}
// TEST_F( TestCategory, batched_scalar_serial_gemm_ct_nt_dcomplex_double ) {
//   typedef ::Test::Gemm::ParamTag<Trans::ConjTranspose,Trans::NoTranspose>
//   param_tag_type; typedef Algo::Gemm::Blocked algo_tag_type;
//   test_batched_gemm<TestDevice,Kokkos::complex<double>,double,param_tag_type,algo_tag_type>();
// }
// TEST_F( TestCategory, batched_scalar_serial_gemm_nt_ct_dcomplex_double ) {
//   typedef ::Test::Gemm::ParamTag<Trans::NoTranspose,Trans::ConjTranspose>
//   param_tag_type; typedef Algo::Gemm::Blocked algo_tag_type;
//   test_batched_gemm<TestDevice,Kokkos::complex<double>,double,param_tag_type,algo_tag_type>();
// }

#endif
