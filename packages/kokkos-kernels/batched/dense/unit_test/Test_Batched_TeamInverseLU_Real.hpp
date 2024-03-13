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
