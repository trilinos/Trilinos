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

#include "KokkosSparse_ccs2crs.hpp"
#include "KokkosSparse_crs2ccs.hpp"
#include "KokkosKernels_TestUtils.hpp"

namespace Test {
template <class CrsType, class IdType, class MapType, class ValsType, class ColsType>
void check_crs_matrix(CrsType crsMat, IdType ccs_row_ids_d, MapType ccs_col_map_d, ValsType ccs_vals_d, ColsType cols) {
  using ordinal_type = typename CrsType::ordinal_type;
  using size_type    = typename CrsType::size_type;

  using ViewTypeRowIds = decltype(ccs_row_ids_d);
  using ViewTypeColMap = decltype(ccs_col_map_d);
  using ViewTypeVals   = decltype(ccs_vals_d);

  // Copy to host
  typename ViewTypeRowIds::HostMirror ccs_row_ids = Kokkos::create_mirror_view(ccs_row_ids_d);
  Kokkos::deep_copy(ccs_row_ids, ccs_row_ids_d);
  typename ViewTypeColMap::HostMirror ccs_col_map = Kokkos::create_mirror_view(ccs_col_map_d);
  Kokkos::deep_copy(ccs_col_map, ccs_col_map_d);
  typename ViewTypeVals::HostMirror ccs_vals = Kokkos::create_mirror_view(ccs_vals_d);
  Kokkos::deep_copy(ccs_vals, ccs_vals_d);

  auto crs_col_ids_d = crsMat.graph.entries;
  auto crs_row_map_d = crsMat.graph.row_map;
  auto crs_vals_d    = crsMat.values;

  using ViewTypeCrsColIds = decltype(crs_col_ids_d);
  using ViewTypeCrsRowMap = decltype(crs_row_map_d);
  using ViewTypeCrsVals   = decltype(crs_vals_d);

  // Copy to host
  typename ViewTypeCrsColIds::HostMirror crs_col_ids = Kokkos::create_mirror_view(crs_col_ids_d);
  Kokkos::deep_copy(crs_col_ids, crs_col_ids_d);
  typename ViewTypeCrsRowMap::HostMirror crs_row_map = Kokkos::create_mirror_view(crs_row_map_d);
  Kokkos::deep_copy(crs_row_map, crs_row_map_d);
  typename ViewTypeCrsVals::HostMirror crs_vals = Kokkos::create_mirror_view(crs_vals_d);
  Kokkos::deep_copy(crs_vals, crs_vals_d);

  Kokkos::fence();

  for (ordinal_type j = 0; j < cols; ++j) {
    auto col_start = ccs_col_map(j);
    auto col_len   = ccs_col_map(j + 1) - col_start;

    for (size_type k = 0; k < col_len; ++k) {
      auto i = col_start + k;

      auto row_start = crs_row_map(ccs_row_ids(i));
      auto row_len   = crs_row_map(ccs_row_ids(i) + 1) - row_start;
      auto row_end   = row_start + row_len;

      if (row_len == 0) continue;

      // Linear search for corresponding element in crs matrix
      auto l = row_start;
      while (l < row_end && crs_col_ids(l) != j) {
        ++l;
      }

      if (l == row_end)
        FAIL() << "crs element at (i: " << ccs_row_ids(i) << ", j: " << j << ") not found!" << std::endl;

      ASSERT_EQ(ccs_vals(i), crs_vals(l)) << "(i: " << ccs_row_ids(i) << ", j: " << j << ")" << std::endl;
    }
  }
}
template <class ScalarType, class LayoutType, class ExeSpaceType>
void doCcs2Crs(size_t m, size_t n, ScalarType min_val, ScalarType max_val, bool fully_sparse = false) {
  RandCsMatrix<ScalarType, LayoutType, ExeSpaceType> ccsMat(n, m, min_val, max_val, fully_sparse);

  auto crsMat = KokkosSparse::ccs2crs(ccsMat.get_dim2(), ccsMat.get_dim1(), ccsMat.get_nnz(), ccsMat.get_vals(),
                                      ccsMat.get_map(), ccsMat.get_ids());

  auto ccs_row_ids_d = ccsMat.get_ids();
  auto ccs_col_map_d = ccsMat.get_map();
  auto ccs_vals_d    = ccsMat.get_vals();
  auto cols          = ccsMat.get_dim1();

  check_crs_matrix(crsMat, ccs_row_ids_d, ccs_col_map_d, ccs_vals_d, cols);
}

template <class LayoutType, class ExeSpaceType>
void doAllScalarsCcs2Crs(size_t m, size_t n, int min, int max) {
  doCcs2Crs<float, LayoutType, ExeSpaceType>(m, n, min, max);
  doCcs2Crs<double, LayoutType, ExeSpaceType>(m, n, min, max);
  doCcs2Crs<Kokkos::complex<float>, LayoutType, ExeSpaceType>(m, n, min, max);
  doCcs2Crs<Kokkos::complex<double>, LayoutType, ExeSpaceType>(m, n, min, max);
}

template <class ExeSpaceType>
void doAllLayoutsCcs2Crs(size_t m, size_t n, int min, int max) {
  doAllScalarsCcs2Crs<Kokkos::LayoutLeft, ExeSpaceType>(m, n, min, max);
  doAllScalarsCcs2Crs<Kokkos::LayoutRight, ExeSpaceType>(m, n, min, max);
}

template <class ExeSpaceType>
void doAllCcs2crs(size_t m, size_t n) {
  int min = 1, max = 10;
  doAllLayoutsCcs2Crs<ExeSpaceType>(m, n, min, max);
}

TEST_F(TestCategory, sparse_ccs2crs) {
  uint64_t ticks = std::chrono::high_resolution_clock::now().time_since_epoch().count() % UINT32_MAX;
  std::srand(ticks);

  // Empty cases
  doCcs2Crs<float, Kokkos::LayoutLeft, TestDevice>(1, 0, 1, 10);
  doCcs2Crs<float, Kokkos::LayoutLeft, TestDevice>(0, 1, 1, 10);

  doCcs2Crs<float, Kokkos::LayoutRight, TestDevice>(1, 0, 1, 10);
  doCcs2Crs<float, Kokkos::LayoutRight, TestDevice>(0, 1, 1, 10);

  doCcs2Crs<float, Kokkos::LayoutLeft, TestDevice>(0, 0, 1, 10);
  doCcs2Crs<float, Kokkos::LayoutRight, TestDevice>(0, 0, 1, 10);

  // Square cases
  for (size_t i = 4; i < 1024; i *= 4) {
    size_t dim = (std::rand() % 511) + 1;
    doAllCcs2crs<TestDevice>(dim, dim);
  }

  // Non-square cases
  for (size_t i = 1; i < 1024; i *= 4) {
    size_t m = (std::rand() % 511) + 1;
    size_t n = (std::rand() % 511) + 1;
    while (n == m) n = (std::rand() % 511) + 1;
    doAllCcs2crs<TestDevice>(m, n);
  }

  // Fully sparse cases
  doCcs2Crs<float, Kokkos::LayoutLeft, TestDevice>(5, 5, 1, 10, true);
  doCcs2Crs<double, Kokkos::LayoutRight, TestDevice>(50, 10, 10, 100, true);

  // Test the convenience wrapper that accepts a ccs matrix
  RandCsMatrix<double, Kokkos::LayoutRight, TestDevice> csMat(2, 2, 10, 10, false);
  auto ccsMatrix =
      crs2ccs(csMat.get_dim1(), csMat.get_dim2(), csMat.get_nnz(), csMat.get_vals(), csMat.get_map(), csMat.get_ids());
  auto crsMatrix = ccs2crs(ccsMatrix);

  auto ccs_row_ids_d = ccsMatrix.graph.entries;
  auto ccs_col_map_d = ccsMatrix.graph.col_map;
  auto ccs_vals_d    = ccsMatrix.values;
  auto cols          = ccsMatrix.numCols();

  check_crs_matrix(crsMatrix, ccs_row_ids_d, ccs_col_map_d, ccs_vals_d, cols);
}
}  // namespace Test