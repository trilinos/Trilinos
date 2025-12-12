// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_solvelu_dcomplex) {
  // printf("Batched serial solveLU - double complex - algorithm type:
  // Unblocked\n");
  test_batched_solvelu<TestDevice, Kokkos::complex<double>, Algo::SolveLU::Unblocked>();
  // printf("Batched serial solveLU - double complex - algorithm type:
  // Blocked\n");
  test_batched_solvelu<TestDevice, Kokkos::complex<double>, Algo::SolveLU::Blocked>();
}
#endif
