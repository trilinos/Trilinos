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

#include "KokkosKernels_TestUtils.hpp"

namespace Test {
template <class ScalarType, class LayoutType, class ExeSpaceType>
void doCsMat(size_t m, size_t n, ScalarType min_val, ScalarType max_val) {
  auto expected_min    = ScalarType(1.0);
  int64_t expected_nnz = 0;
  RandCsMatrix<ScalarType, LayoutType, ExeSpaceType> cm(m, n, min_val, max_val);

  for (int64_t i = 0; i < cm.get_nnz(); ++i)
    ASSERT_GE(cm(i), expected_min) << cm.info;

  auto map_d = cm.get_map();
  auto map   = Kokkos::create_mirror_view(map_d);
  Kokkos::deep_copy(map, map_d);

  for (int64_t j = 0; j < cm.get_dim1(); ++j) {
    int64_t row_len = j < static_cast<int64_t>(m) ? (map(j + 1) - map(j)) : 0;
    for (int64_t i = 0; i < row_len; ++i) {
      int64_t row_start = j < static_cast<int64_t>(m) ? map(j) : 0;
      ASSERT_FLOAT_EQ(cm(row_start + i), cm(expected_nnz + i)) << cm.info;
    }
    expected_nnz += row_len;
  }
  ASSERT_EQ(cm.get_nnz(), expected_nnz) << cm.info;

  // No need to check data here. Kokkos unit-tests deep_copy.
  auto vals = cm.get_vals();
  ASSERT_EQ(vals.extent(0), cm.get_nnz() + 1) << cm.info;

  auto row_ids = cm.get_ids();
  ASSERT_EQ(row_ids.extent(0), cm.get_dim1() * cm.get_dim2() + 1) << cm.info;

  auto col_map = cm.get_map();
  ASSERT_EQ(col_map.extent(0), cm.get_dim1() + 1);
}

template <class ExeSpaceType>
void doAllCsMat(size_t m, size_t n) {
  int min = 1, max = 10;

  // Verify that CsMax is constructed properly.
  doCsMat<float, Kokkos::LayoutLeft, ExeSpaceType>(m, n, min, max);
  doCsMat<float, Kokkos::LayoutRight, ExeSpaceType>(m, n, min, max);

  doCsMat<double, Kokkos::LayoutLeft, ExeSpaceType>(m, n, min, max);
  doCsMat<double, Kokkos::LayoutRight, ExeSpaceType>(m, n, min, max);

  // Verify that CsMat can be instantiated with complex types.
  RandCsMatrix<Kokkos::complex<float>, Kokkos::LayoutLeft, ExeSpaceType> cmcf(
      m, n, min, max);
  RandCsMatrix<Kokkos::complex<double>, Kokkos::LayoutRight, ExeSpaceType> cmcd(
      m, n, min, max);
}

// Test randomly generated Cs matrices
TEST_F(TestCategory, sparse_randcsmat) {
  // Square cases
  for (int dim = 1; dim < 1024; dim *= 4) doAllCsMat<TestExecSpace>(dim, dim);

  // Non-square cases
  for (int dim = 1; dim < 1024; dim *= 4) {
    doAllCsMat<TestExecSpace>(dim * 3, dim);
    doAllCsMat<TestExecSpace>(dim, dim * 3);
  }
}
}  // namespace Test