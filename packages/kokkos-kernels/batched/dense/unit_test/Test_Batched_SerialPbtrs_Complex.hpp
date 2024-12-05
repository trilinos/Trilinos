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
TEST_F(TestCategory, test_batched_pbtrs_l_fcomplex) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<Uplo::Lower>;

  test_batched_pbtrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pbtrs_u_fcomplex) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<Uplo::Upper>;

  test_batched_pbtrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_pbtrs_l_dcomplex) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<Uplo::Lower>;

  test_batched_pbtrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pbtrs_u_dcomplex) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<Uplo::Upper>;

  test_batched_pbtrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
#endif
