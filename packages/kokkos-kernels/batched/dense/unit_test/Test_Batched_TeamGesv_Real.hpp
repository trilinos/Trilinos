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
TEST_F(TestCategory, batched_scalar_team_gesv_static_pivoting_float) {
  test_batched_team_gesv<TestDevice, float, KokkosBatched::Gesv::StaticPivoting>();
}
TEST_F(TestCategory, batched_scalar_team_gesv_no_pivoting_float) {
  test_batched_team_gesv<TestDevice, float, KokkosBatched::Gesv::NoPivoting>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_team_gesv_static_pivoting_double) {
  test_batched_team_gesv<TestDevice, double, KokkosBatched::Gesv::StaticPivoting>();
}
TEST_F(TestCategory, batched_scalar_team_gesv_no_pivoting_double) {
  test_batched_team_gesv<TestDevice, double, KokkosBatched::Gesv::NoPivoting>();
}
#endif
