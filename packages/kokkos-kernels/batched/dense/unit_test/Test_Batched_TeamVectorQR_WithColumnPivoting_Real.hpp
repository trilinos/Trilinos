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
TEST_F(TestCategory, batched_scalar_teamvector_qr_with_columnpivoting_float) {
  typedef Algo::QR::Unblocked algo_tag_type;
  test_batched_qr_with_columnpivoting<TestDevice, float, int, algo_tag_type>();
}
#endif

// FIXME_SYCL timeout
#ifndef KOKKOS_ENABLE_SYCL
#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_teamvector_qr_with_columnpivoting_double) {
  typedef Algo::QR::Unblocked algo_tag_type;
  test_batched_qr_with_columnpivoting<TestDevice, double, int, algo_tag_type>();
}
#endif
#endif
