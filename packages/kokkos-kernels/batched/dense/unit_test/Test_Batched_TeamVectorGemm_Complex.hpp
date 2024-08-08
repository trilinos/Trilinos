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
TEST_F(TestCategory, batched_scalar_team_vector_gemm_nt_nt_scomplex_scomplex) {
  typedef ::Test::TeamVectorGemm::ParamTag<Trans::NoTranspose, Trans::NoTranspose> param_tag_type;

  // test_batched_teamvectorgemm<TestDevice,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                              Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_vector_gemm_t_nt_scomplex_scomplex) {
  typedef ::Test::TeamVectorGemm::ParamTag<Trans::Transpose, Trans::NoTranspose> param_tag_type;

  // test_batched_teamvectorgemm<TestDevice,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                              Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_vector_gemm_nt_t_scomplex_scomplex) {
  typedef ::Test::TeamVectorGemm::ParamTag<Trans::NoTranspose, Trans::Transpose> param_tag_type;

  // test_batched_teamvectorgemm<TestDevice,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                              Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_vector_gemm_t_t_scomplex_scomplex) {
  typedef ::Test::TeamVectorGemm::ParamTag<Trans::Transpose, Trans::Transpose> param_tag_type;

  // test_batched_teamvectorgemm<TestDevice,Kokkos::complex<float>,Kokkos::complex<float>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>, param_tag_type,
                              Algo::Gemm::Unblocked>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, batched_scalar_team_vector_gemm_nt_nt_dcomplex_dcomplex) {
  typedef ::Test::TeamVectorGemm::ParamTag<Trans::NoTranspose, Trans::NoTranspose> param_tag_type;

  // test_batched_teamvectorgemm<TestDevice,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                              Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_vector_gemm_t_nt_dcomplex_dcomplex) {
  typedef ::Test::TeamVectorGemm::ParamTag<Trans::Transpose, Trans::NoTranspose> param_tag_type;

  // test_batched_teamvectorgemm<TestDevice,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                              Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_vector_gemm_nt_t_dcomplex_dcomplex) {
  typedef ::Test::TeamVectorGemm::ParamTag<Trans::NoTranspose, Trans::Transpose> param_tag_type;

  // test_batched_teamvectorgemm<TestDevice,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                              Algo::Gemm::Unblocked>();
}
TEST_F(TestCategory, batched_scalar_team_vector_gemm_t_t_dcomplex_dcomplex) {
  typedef ::Test::TeamVectorGemm::ParamTag<Trans::Transpose, Trans::Transpose> param_tag_type;

  // test_batched_teamvectorgemm<TestDevice,Kokkos::complex<double>,Kokkos::complex<double>,param_tag_type,Algo::Gemm::Blocked>();
  test_batched_teamvectorgemm<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>, param_tag_type,
                              Algo::Gemm::Unblocked>();
}
#endif
