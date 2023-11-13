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
TEST_F(TestCategory, batched_scalar_teamvector_solve_utv_float) {
  typedef Algo::UTV::Unblocked algo_tag_type;
  test_batched_solve_utv<TestExecSpace, float, int, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// FIXME_SYCL
#ifndef KOKKOS_ENABLE_SYCL
TEST_F(TestCategory, batched_scalar_teamvector_solve_utv_double) {
  typedef Algo::UTV::Unblocked algo_tag_type;
  test_batched_solve_utv<TestExecSpace, double, int, algo_tag_type>();
}
#endif
#endif
