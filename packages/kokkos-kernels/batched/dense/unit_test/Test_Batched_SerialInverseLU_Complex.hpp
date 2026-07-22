// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_inverselu_dcomplex) {
  // printf("Batched serial inverse LU - double complex - algorithm type:
  // Unblocked\n");
  test_batched_inverselu<TestDevice, Kokkos::complex<double>, Algo::InverseLU::Unblocked>();
  // printf("Batched serial inverse LU - double complex - algorithm type:
  // Blocked\n");
  test_batched_inverselu<TestDevice, Kokkos::complex<double>, Algo::InverseLU::Blocked>();
}
#endif
