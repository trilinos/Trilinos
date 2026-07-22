// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, batched_scalar_team_lu_dcomplex) {
  typedef Algo::LU::Blocked algo_tag_type;
  test_batched_lu<TestDevice, Kokkos::complex<double>, algo_tag_type>();
}
#endif
