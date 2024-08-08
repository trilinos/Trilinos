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

#include "KokkosSparse_coo2crs.hpp"
#include "KokkosSparse_crs2coo.hpp"
#include "KokkosKernels_TestUtils.hpp"

namespace Test {
template <class CrsType, class RowType, class ColType, class DataType>
void check_coo_matrix(CrsType crsMatRef, RowType row, ColType col, DataType data) {
  // Copy coo to host
  typename RowType::HostMirror row_h = Kokkos::create_mirror_view(row);
  Kokkos::deep_copy(row_h, row);
  typename ColType::HostMirror col_h = Kokkos::create_mirror_view(col);
  Kokkos::deep_copy(col_h, col);
  typename DataType::HostMirror data_h = Kokkos::create_mirror_view(data);
  Kokkos::deep_copy(data_h, data);

  // printf("coo in:\n");
  // for (unsigned i = 0; i < data_h.extent(0); i++)
  // printf("(%lld, %lld, %g)\n", row_h(i), col_h(i), data_h(i));

  auto crs_col_ids_ref_d = crsMatRef.graph.entries;
  auto crs_row_map_ref_d = crsMatRef.graph.row_map;
  auto crs_vals_ref_d    = crsMatRef.values;

  using ViewTypeCrsColIdsRef = decltype(crs_col_ids_ref_d);
  using ViewTypeCrsRowMapRef = decltype(crs_row_map_ref_d);
  using ViewTypeCrsValsRef   = decltype(crs_vals_ref_d);

  // Copy crs to host
  typename ViewTypeCrsColIdsRef::HostMirror crs_col_ids_ref = Kokkos::create_mirror_view(crs_col_ids_ref_d);
  Kokkos::deep_copy(crs_col_ids_ref, crs_col_ids_ref_d);
  typename ViewTypeCrsRowMapRef::HostMirror crs_row_map_ref = Kokkos::create_mirror_view(crs_row_map_ref_d);
  Kokkos::deep_copy(crs_row_map_ref, crs_row_map_ref_d);
  typename ViewTypeCrsValsRef::HostMirror crs_vals_ref = Kokkos::create_mirror_view(crs_vals_ref_d);
  Kokkos::deep_copy(crs_vals_ref, crs_vals_ref_d);

  Kokkos::fence();

  ASSERT_EQ(crsMatRef.nnz(), row.extent(0));
  ASSERT_EQ(crsMatRef.nnz(), col.extent(0));
  ASSERT_EQ(crsMatRef.nnz(), data.extent(0));

  for (decltype(row.extent(0)) idx = 0; idx < row.extent(0); ++idx) {
    auto row_id = row_h(idx);
    auto col_id = col_h(idx);
    auto val    = data_h(idx);
    std::string fail_msg =
        "idx - " + std::to_string(idx) + " row: " + std::to_string(row_id) + ", col: " + std::to_string(col_id);

    auto row_start_ref = crs_row_map_ref(row_id);
    auto row_stop_ref  = crs_row_map_ref(row_id + 1);

    auto crs_idx = row_start_ref;
    for (; crs_idx < row_stop_ref; crs_idx++) {
      if (crs_col_ids_ref(crs_idx) == col_id) {
        // crs2coo does a direct copy, no need for an epsilon.
        if (crs_vals_ref(crs_idx) == val) break;
      }
    }
    if (crs_idx == row_stop_ref) FAIL() << fail_msg << " not found in crsMatRef!";
  }
}

template <class ScalarType, class LayoutType, class ExeSpaceType>
void doCrs2Coo(size_t m, size_t n, ScalarType min_val, ScalarType max_val) {
  using RandCrsMatType = RandCsMatrix<ScalarType, LayoutType, ExeSpaceType>;
  RandCrsMatType crsMat(m, n, min_val, max_val, m == 0 || n == 0);

  using CrsOT   = typename RandCrsMatType::IdViewTypeD::value_type;
  using CrsType = typename KokkosSparse::CrsMatrix<ScalarType, CrsOT, ExeSpaceType>;
  auto map      = crsMat.get_map();
  auto ids      = crsMat.get_ids();
  CrsType crsMatrix("doCrs2Coo", crsMat.get_dim1(), crsMat.get_dim2(), crsMat.get_nnz(), crsMat.get_vals(), map, ids);

  auto cooMat = KokkosSparse::crs2coo(crsMatrix);
  check_coo_matrix(crsMatrix, cooMat.row(), cooMat.col(), cooMat.data());
}

template <class LayoutType, class ExeSpaceType>
void doAllScalarsCrs2Coo(size_t m, size_t n, int min, int max) {
  doCrs2Coo<float, LayoutType, ExeSpaceType>(m, n, min, max);
  doCrs2Coo<double, LayoutType, ExeSpaceType>(m, n, min, max);
  doCrs2Coo<Kokkos::complex<float>, LayoutType, ExeSpaceType>(m, n, min, max);
  doCrs2Coo<Kokkos::complex<double>, LayoutType, ExeSpaceType>(m, n, min, max);
}

template <class ExeSpaceType>
void doAllLayoutsCrs2Coo(size_t m, size_t n, int min, int max) {
  doAllScalarsCrs2Coo<Kokkos::LayoutLeft, ExeSpaceType>(m, n, min, max);
  doAllScalarsCrs2Coo<Kokkos::LayoutRight, ExeSpaceType>(m, n, min, max);
}

template <class ExeSpaceType>
void doAllCrs2Coo(size_t m, size_t n) {
  int min = 1, max = 10;
  doAllLayoutsCrs2Coo<ExeSpaceType>(m, n, min, max);
}

TEST_F(TestCategory, sparse_crs2coo) {
  uint64_t ticks = std::chrono::high_resolution_clock::now().time_since_epoch().count() % UINT32_MAX;
  std::srand(ticks);

  // Square cases
  for (size_t i = 1; i < 256; i *= 4) {
    size_t dim = (std::rand() % 511) + 1;
    doAllCrs2Coo<TestDevice>(dim, dim);
  }

  // Non-square cases
  for (size_t i = 1; i < 256; i *= 4) {
    size_t m = (std::rand() % 511) + 1;
    size_t n = (std::rand() % 511) + 1;
    while (n == m) n = (std::rand() % 511) + 1;
    doAllCrs2Coo<TestDevice>(m, n);
  }
}
}  // namespace Test