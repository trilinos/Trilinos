// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_solvelu_float) {
  // printf("Batched serial solveLU - float - algorithm type: Unblocked\n");
  test_batched_solvelu<TestDevice, float, Algo::SolveLU::Unblocked>();
  // printf("Batched serial solveLU - float - algorithm type: Blocked\n");
  test_batched_solvelu<TestDevice, float, Algo::SolveLU::Blocked>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_solvelu_double) {
  // printf("Batched serial solveLU - double - algorithm type: Unblocked\n");
  test_batched_solvelu<TestDevice, double, Algo::SolveLU::Unblocked>();
  // printf("Batched serial solveLU - double - algorithm type: Blocked\n");
  test_batched_solvelu<TestDevice, double, Algo::SolveLU::Blocked>();
}
#endif
