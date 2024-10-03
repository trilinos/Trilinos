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

#include "KokkosSparse_csc2csr.hpp"
#include "KokkosKernels_TestUtils.hpp"

namespace Test {
template <class ScalarType, class LayoutType, class ExeSpaceType>
void doCsc2Csr(size_t m, size_t n, ScalarType min_val, ScalarType max_val, bool fully_sparse = false) {
  RandCsMatrix<ScalarType, LayoutType, ExeSpaceType> cscMat(n, m, min_val, max_val, fully_sparse);

  auto csrMat = KokkosSparse::csc2csr(cscMat.get_dim2(), cscMat.get_dim1(), cscMat.get_nnz(), cscMat.get_vals(),
                                      cscMat.get_map(), cscMat.get_ids());

  auto csc_row_ids_d = cscMat.get_ids();
  auto csc_col_map_d = cscMat.get_map();
  auto csc_vals_d    = cscMat.get_vals();

  using ViewTypeRowIds = decltype(csc_row_ids_d);
  using ViewTypeColMap = decltype(csc_col_map_d);
  using ViewTypeVals   = decltype(csc_vals_d);

  // Copy to host
  typename ViewTypeRowIds::HostMirror csc_row_ids = Kokkos::create_mirror_view(csc_row_ids_d);
  Kokkos::deep_copy(csc_row_ids, csc_row_ids_d);
  typename ViewTypeColMap::HostMirror csc_col_map = Kokkos::create_mirror_view(csc_col_map_d);
  Kokkos::deep_copy(csc_col_map, csc_col_map_d);
  typename ViewTypeVals::HostMirror csc_vals = Kokkos::create_mirror_view(csc_vals_d);
  Kokkos::deep_copy(csc_vals, csc_vals_d);

  auto csr_col_ids_d = csrMat.graph.entries;
  auto csr_row_map_d = csrMat.graph.row_map;
  auto csr_vals_d    = csrMat.values;

  using ViewTypeCsrColIds = decltype(csr_col_ids_d);
  using ViewTypeCsrRowMap = decltype(csr_row_map_d);
  using ViewTypeCsrVals   = decltype(csr_vals_d);

  // Copy to host
  typename ViewTypeCsrColIds::HostMirror csr_col_ids = Kokkos::create_mirror_view(csr_col_ids_d);
  Kokkos::deep_copy(csr_col_ids, csr_col_ids_d);
  typename ViewTypeCsrRowMap::HostMirror csr_row_map = Kokkos::create_mirror_view(csr_row_map_d);
  Kokkos::deep_copy(csr_row_map, csr_row_map_d);
  typename ViewTypeCsrVals::HostMirror csr_vals = Kokkos::create_mirror_view(csr_vals_d);
  Kokkos::deep_copy(csr_vals, csr_vals_d);

  Kokkos::fence();

  for (int j = 0; j < cscMat.get_dim1(); ++j) {
    auto col_start = csc_col_map(j);
    auto col_len   = csc_col_map(j + 1) - col_start;

    for (int k = 0; k < col_len; ++k) {
      auto i = col_start + k;

      auto row_start = csr_row_map(csc_row_ids(i));
      auto row_len   = csr_row_map(csc_row_ids(i) + 1) - row_start;
      auto row_end   = row_start + row_len;

      if (row_len == 0) continue;

      // Linear search for corresponding element in csr matrix
      int l = row_start;
      while (l < row_end && csr_col_ids(l) != j) {
        ++l;
      }

      if (l == row_end)
        FAIL() << "csr element at (i: " << csc_row_ids(i) << ", j: " << j << ") not found!" << std::endl;

      ASSERT_EQ(csc_vals(i), csr_vals(l)) << "(i: " << csc_row_ids(i) << ", j: " << j << ")" << std::endl;
    }
  }
}

template <class LayoutType, class ExeSpaceType>
void doAllScalarsCsc2Csr(size_t m, size_t n, int min, int max) {
  doCsc2Csr<float, LayoutType, ExeSpaceType>(m, n, min, max);
  doCsc2Csr<double, LayoutType, ExeSpaceType>(m, n, min, max);
  doCsc2Csr<Kokkos::complex<float>, LayoutType, ExeSpaceType>(m, n, min, max);
  doCsc2Csr<Kokkos::complex<double>, LayoutType, ExeSpaceType>(m, n, min, max);
}

template <class ExeSpaceType>
void doAllLayoutsCsc2Csr(size_t m, size_t n, int min, int max) {
  doAllScalarsCsc2Csr<Kokkos::LayoutLeft, ExeSpaceType>(m, n, min, max);
  doAllScalarsCsc2Csr<Kokkos::LayoutRight, ExeSpaceType>(m, n, min, max);
}

template <class ExeSpaceType>
void doAllCsc2csr(size_t m, size_t n) {
  int min = 1, max = 10;
  doAllLayoutsCsc2Csr<ExeSpaceType>(m, n, min, max);
}

TEST_F(TestCategory, sparse_csc2csr) {
  uint64_t ticks = std::chrono::high_resolution_clock::now().time_since_epoch().count() % UINT32_MAX;
  std::srand(ticks);

  // Empty cases
  doCsc2Csr<float, Kokkos::LayoutLeft, TestDevice>(1, 0, 1, 10);
  doCsc2Csr<float, Kokkos::LayoutLeft, TestDevice>(0, 1, 1, 10);

  doCsc2Csr<float, Kokkos::LayoutRight, TestDevice>(1, 0, 1, 10);
  doCsc2Csr<float, Kokkos::LayoutRight, TestDevice>(0, 1, 1, 10);

  doCsc2Csr<float, Kokkos::LayoutLeft, TestDevice>(0, 0, 1, 10);
  doCsc2Csr<float, Kokkos::LayoutRight, TestDevice>(0, 0, 1, 10);

  // Square cases
  for (size_t i = 4; i < 1024; i *= 4) {
    size_t dim = (std::rand() % 511) + 1;
    doAllCsc2csr<TestDevice>(dim, dim);
  }

  // Non-square cases
  for (size_t i = 1; i < 1024; i *= 4) {
    size_t m = (std::rand() % 511) + 1;
    size_t n = (std::rand() % 511) + 1;
    while (n == m) n = (std::rand() % 511) + 1;
    doAllCsc2csr<TestDevice>(m, n);
  }

  // Fully sparse cases
  doCsc2Csr<float, Kokkos::LayoutLeft, TestDevice>(5, 5, 1, 10, true);
  doCsc2Csr<double, Kokkos::LayoutRight, TestDevice>(50, 10, 10, 100, true);
}
}  // namespace Test