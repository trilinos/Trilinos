// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_team_inverselu_float) {
  // printf("Batched team inverse LU - float - algorithm type: Unblocked\n");
  test_batched_inverselu<TestDevice, float, Algo::InverseLU::Unblocked>();
  // printf("Batched team inverse LU - float - algorithm type: Blocked\n");
  test_batched_inverselu<TestDevice, float, Algo::InverseLU::Blocked>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_team_inverselu_double) {
  // printf("Batched team inverse LU - double - algorithm type: Unblocked\n");
  test_batched_inverselu<TestDevice, double, Algo::InverseLU::Unblocked>();
  // printf("Batched team inverse LU - double - algorithm type: Blocked\n");
  test_batched_inverselu<TestDevice, double, Algo::InverseLU::Blocked>();
}
#endif
