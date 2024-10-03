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

#include "KokkosSparse_crs2ccs.hpp"
#include "KokkosSparse_ccs2crs.hpp"
#include "KokkosKernels_TestUtils.hpp"

namespace Test {
template <class CcsType, class IdType, class MapType, class ValsType, class ColsType>
void check_ccs_matrix(CcsType ccsMat, IdType crs_col_ids_d, MapType crs_row_map_d, ValsType crs_vals_d, ColsType cols) {
  using ordinal_type = typename CcsType::ordinal_type;
  using size_type    = typename CcsType::size_type;

  using ViewTypeRowIds = decltype(crs_col_ids_d);
  using ViewTypeColMap = decltype(crs_row_map_d);
  using ViewTypeVals   = decltype(crs_vals_d);

  // Copy to host
  typename ViewTypeRowIds::HostMirror crs_col_ids = Kokkos::create_mirror_view(crs_col_ids_d);
  Kokkos::deep_copy(crs_col_ids, crs_col_ids_d);
  typename ViewTypeColMap::HostMirror crs_row_map = Kokkos::create_mirror_view(crs_row_map_d);
  Kokkos::deep_copy(crs_row_map, crs_row_map_d);
  typename ViewTypeVals::HostMirror crs_vals = Kokkos::create_mirror_view(crs_vals_d);
  Kokkos::deep_copy(crs_vals, crs_vals_d);

  auto ccs_row_ids_d = ccsMat.graph.entries;
  auto ccs_col_map_d = ccsMat.graph.col_map;
  auto ccs_vals_d    = ccsMat.values;

  using ViewTypeCrsColIds = decltype(ccs_row_ids_d);
  using ViewTypeCrsRowMap = decltype(ccs_col_map_d);
  using ViewTypeCrsVals   = decltype(ccs_vals_d);

  // Copy to host
  typename ViewTypeCrsColIds::HostMirror ccs_row_ids = Kokkos::create_mirror_view(ccs_row_ids_d);
  Kokkos::deep_copy(ccs_row_ids, ccs_row_ids_d);
  typename ViewTypeCrsRowMap::HostMirror ccs_col_map = Kokkos::create_mirror_view(ccs_col_map_d);
  Kokkos::deep_copy(ccs_col_map, ccs_col_map_d);
  typename ViewTypeCrsVals::HostMirror ccs_vals = Kokkos::create_mirror_view(ccs_vals_d);
  Kokkos::deep_copy(ccs_vals, ccs_vals_d);

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
        FAIL() << "ccs element at (i: " << ccs_row_ids(i) << ", j: " << j << ") not found!" << std::endl;

      ASSERT_EQ(ccs_vals(i), crs_vals(l)) << "(i: " << ccs_row_ids(i) << ", j: " << j << ")" << std::endl;
    }
  }
}

template <class ScalarType, class LayoutType, class ExeSpaceType>
void doCrs2Ccs(size_t m, size_t n, ScalarType min_val, ScalarType max_val, bool fully_sparse = false) {
  RandCsMatrix<ScalarType, LayoutType, ExeSpaceType> crsMat(m, n, min_val, max_val, fully_sparse);

  auto ccsMat = KokkosSparse::crs2ccs(crsMat.get_dim1(), crsMat.get_dim2(), crsMat.get_nnz(), crsMat.get_vals(),
                                      crsMat.get_map(), crsMat.get_ids());

  auto crs_col_ids_d = crsMat.get_ids();
  auto crs_row_map_d = crsMat.get_map();
  auto crs_vals_d    = crsMat.get_vals();
  auto cols          = crsMat.get_dim2();
  check_ccs_matrix(ccsMat, crs_col_ids_d, crs_row_map_d, crs_vals_d, cols);
}

template <class LayoutType, class ExeSpaceType>
void doAllScalarsCrs2Ccs(size_t m, size_t n, int min, int max) {
  doCrs2Ccs<float, LayoutType, ExeSpaceType>(m, n, min, max);
  doCrs2Ccs<double, LayoutType, ExeSpaceType>(m, n, min, max);
  doCrs2Ccs<Kokkos::complex<float>, LayoutType, ExeSpaceType>(m, n, min, max);
  doCrs2Ccs<Kokkos::complex<double>, LayoutType, ExeSpaceType>(m, n, min, max);
}

template <class ExeSpaceType>
void doAllLayoutsCrs2Ccs(size_t m, size_t n, int min, int max) {
  doAllScalarsCrs2Ccs<Kokkos::LayoutLeft, ExeSpaceType>(m, n, min, max);
  doAllScalarsCrs2Ccs<Kokkos::LayoutRight, ExeSpaceType>(m, n, min, max);
}

template <class ExeSpaceType>
void doAllCrs2Ccs(size_t m, size_t n) {
  int min = 1, max = 10;
  doAllLayoutsCrs2Ccs<ExeSpaceType>(m, n, min, max);
}

TEST_F(TestCategory, sparse_crs2ccs) {
  uint64_t ticks = std::chrono::high_resolution_clock::now().time_since_epoch().count() % UINT32_MAX;
  std::srand(ticks);

  // Empty cases
  doCrs2Ccs<float, Kokkos::LayoutLeft, TestDevice>(1, 0, 1, 10);
  doCrs2Ccs<float, Kokkos::LayoutLeft, TestDevice>(0, 1, 1, 10);

  doCrs2Ccs<float, Kokkos::LayoutRight, TestDevice>(1, 0, 1, 10);
  doCrs2Ccs<float, Kokkos::LayoutRight, TestDevice>(0, 1, 1, 10);

  doCrs2Ccs<float, Kokkos::LayoutLeft, TestDevice>(0, 0, 1, 10);
  doCrs2Ccs<float, Kokkos::LayoutRight, TestDevice>(0, 0, 1, 10);

  // Square cases
  for (size_t i = 4; i < 1024; i *= 4) {
    size_t dim = (std::rand() % 511) + 1;
    doAllCrs2Ccs<TestDevice>(dim, dim);
  }

  // Non-square cases
  for (size_t i = 1; i < 1024; i *= 4) {
    size_t m = (std::rand() % 511) + 1;
    size_t n = (std::rand() % 511) + 1;
    while (n == m) n = (std::rand() % 511) + 1;
    doAllCrs2Ccs<TestDevice>(m, n);
  }

  // Fully sparse cases
  doCrs2Ccs<float, Kokkos::LayoutLeft, TestDevice>(5, 5, 1, 10, true);
  doCrs2Ccs<double, Kokkos::LayoutRight, TestDevice>(50, 10, 10, 100, true);

  // Test the convenience wrapper that accepts a crs matrix
  RandCsMatrix<float, Kokkos::LayoutLeft, TestDevice> csMat(2, 2, 10, 10, false);
  auto crsMatrix =
      ccs2crs(csMat.get_dim2(), csMat.get_dim1(), csMat.get_nnz(), csMat.get_vals(), csMat.get_map(), csMat.get_ids());
  auto ccsMatrix = crs2ccs(crsMatrix);

  auto crs_col_ids_d = crsMatrix.graph.entries;
  auto crs_row_map_d = crsMatrix.graph.row_map;
  auto crs_vals_d    = crsMatrix.values;
  auto cols          = crsMatrix.numCols();
  check_ccs_matrix(ccsMatrix, crs_col_ids_d, crs_row_map_d, crs_vals_d, cols);
}
}  // namespace Test