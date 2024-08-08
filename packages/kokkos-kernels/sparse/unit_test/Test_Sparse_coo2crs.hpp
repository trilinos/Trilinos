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
CrsType vanilla_coo2crs(size_t m, size_t n, RowType row, ColType col, DataType data) {
  using RowIndexType = typename RowType::value_type;
  using ColIndexType = typename ColType::value_type;
  using ValueType    = typename DataType::value_type;
  std::unordered_map<RowIndexType, std::unordered_map<ColIndexType, ValueType> *> umap;
  int nnz = 0;

  for (uint64_t i = 0; i < data.extent(0); i++) {
    auto r = row(i);
    auto c = col(i);
    auto d = data(i);

    if (r >= 0 && c >= 0) {
      if (umap.find(r) != umap.end()) {  // exists
        auto my_row = umap.at(r);
        if (my_row->find(c) != my_row->end())
          my_row->at(c) += d;
        else {
          my_row->insert(std::make_pair(c, d));
          nnz++;
        }
      } else {  // create a new row.
        auto new_row = new std::unordered_map<ColIndexType, ValueType>();
        umap.insert(std::make_pair(r, new_row));
        new_row->insert(std::make_pair(c, d));
        nnz++;
      }
    }
  }

  typename CrsType::row_map_type::non_const_type row_map("vanilla_row_map", m + 1);
  typename CrsType::values_type values("vanilla_values", nnz);
  typename CrsType::staticcrsgraph_type::entries_type col_ids("vanilla_col_ids", nnz);

  typename CrsType::row_map_type::non_const_type::HostMirror row_map_h      = Kokkos::create_mirror_view(row_map);
  typename CrsType::values_type::HostMirror values_h                        = Kokkos::create_mirror_view(values);
  typename CrsType::staticcrsgraph_type::entries_type::HostMirror col_ids_h = Kokkos::create_mirror_view(col_ids);

  int row_len = 0;
  for (uint64_t i = 0; i < m; i++) {
    if (umap.find(i) != umap.end()) row_len += umap.at(i)->size();
    row_map_h(i + 1) = row_len;
  }

  for (uint64_t i = 0; i < m; i++) {
    if (umap.find(i) == umap.end())  // Fully sparse row
      continue;

    auto row_start = row_map_h(i);
    auto row_end   = row_map_h(i + 1);
    auto my_row    = umap.at(i);
    auto iter      = my_row->begin();
    for (auto j = row_start; j < row_end; j++, iter++) {
      col_ids_h(j) = iter->first;
      values_h(j)  = iter->second;
    }
    delete my_row;
  }

  Kokkos::deep_copy(row_map, row_map_h);
  Kokkos::deep_copy(col_ids, col_ids_h);
  Kokkos::deep_copy(values, values_h);

  return CrsType("vanilla_coo2csr", m, n, nnz, values, row_map, col_ids);
}

template <class CrsType, class RowType, class ColType, class DataType>
void check_crs_matrix(CrsType crsMat, RowType row, ColType col, DataType data,
                      std::string failure_info = "no failure information!") {
  using value_type = typename DataType::value_type;
  using ats        = Kokkos::ArithTraits<value_type>;

  // Copy coo to host
  typename RowType::HostMirror row_h = Kokkos::create_mirror_view(row);
  Kokkos::deep_copy(row_h, row);
  typename ColType::HostMirror col_h = Kokkos::create_mirror_view(col);
  Kokkos::deep_copy(col_h, col);
  typename DataType::HostMirror data_h = Kokkos::create_mirror_view(data);
  Kokkos::deep_copy(data_h, data);

  auto crsMatRef =
      vanilla_coo2crs<CrsType, typename RowType::HostMirror, typename ColType::HostMirror,
                      typename DataType::HostMirror>(crsMat.numRows(), crsMat.numCols(), row_h, col_h, data_h);

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

  auto crs_col_ids_d = crsMat.graph.entries;
  auto crs_row_map_d = crsMat.graph.row_map;
  auto crs_vals_d    = crsMat.values;

  using ViewTypeCrsColIds = decltype(crs_col_ids_d);
  using ViewTypeCrsRowMap = decltype(crs_row_map_d);
  using ViewTypeCrsVals   = decltype(crs_vals_d);

  // Copy crs to host
  typename ViewTypeCrsColIds::HostMirror crs_col_ids = Kokkos::create_mirror_view(crs_col_ids_d);
  Kokkos::deep_copy(crs_col_ids, crs_col_ids_d);
  typename ViewTypeCrsRowMap::HostMirror crs_row_map = Kokkos::create_mirror_view(crs_row_map_d);
  Kokkos::deep_copy(crs_row_map, crs_row_map_d);
  typename ViewTypeCrsVals::HostMirror crs_vals = Kokkos::create_mirror_view(crs_vals_d);
  Kokkos::deep_copy(crs_vals, crs_vals_d);

  Kokkos::fence();

  ASSERT_EQ(crsMatRef.nnz(), crsMat.nnz()) << failure_info;

  for (int i = 0; i < crsMatRef.numRows(); i++) {
    ASSERT_EQ(crs_row_map_ref(i), crs_row_map(i))
        << "crs_row_map_ref(" << i << " = " << crs_row_map_ref(i) << " != "
        << "crs_row_map(" << i << " = " << crs_row_map(i) << " -- " << failure_info;
  }

  for (int i = 0; i < crsMatRef.numRows(); ++i) {
    auto row_start_ref = crs_row_map_ref(i);
    auto row_stop_ref  = crs_row_map_ref(i + 1);
    auto row_len_ref   = row_stop_ref - row_start_ref;

    auto row_start = crs_row_map(i);
    auto row_len   = crs_row_map(i + 1) - row_start;

    ASSERT_EQ(row_start_ref, row_start);
    ASSERT_EQ(row_len_ref, row_len);

    for (auto j = row_start_ref; j < row_stop_ref; ++j) {
      // Look for the corresponding col_id
      auto col_id_ref = crs_col_ids_ref(j);
      std::string fail_msg =
          "row: " + std::to_string(i) + ", crs_col_ids_ref(" + std::to_string(j) + ") = " + std::to_string(col_id_ref);

      auto k = row_start_ref;
      for (; k < row_stop_ref; ++k) {
        if (crs_col_ids(k) == col_id_ref) break;
      }
      if (k == row_stop_ref) FAIL() << fail_msg << " not found in crs_col_ids!" << failure_info;

      // NOTE: ASSERT_EQ doesn't work -- values may be summed in different
      // orders We sum at most m x n values.
      auto eps = crsMatRef.numCols() * crsMatRef.numRows() * 10e1 * ats::epsilon();
      EXPECT_NEAR_KK(crs_vals_ref(j), crs_vals(k), eps, fail_msg + " mismatched values!" + failure_info);
    }
  }
}

template <class ScalarType, class LayoutType, class Device>
void doCoo2Crs(size_t m, size_t n, ScalarType min_val, ScalarType max_val) {
  RandCooMat<ScalarType, LayoutType, Device> cooMat(m, n, m * n, min_val, max_val);
  auto randRow  = cooMat.get_row();
  auto randCol  = cooMat.get_col();
  auto randData = cooMat.get_data();

  std::string failure_info = "\nBegin arguments for above failure...\n" + cooMat.info +
                             "scalar: " + std::string(typeid(ScalarType).name()) + "\n" +
                             "layout: " + std::string(typeid(LayoutType).name()) + "\n" + "m: " + std::to_string(m) +
                             ", n: " + std::to_string(n) + "\n...end arguments for above failure.\n";

  auto crsMat = KokkosSparse::coo2crs(m, n, randRow, randCol, randData);
  check_crs_matrix(crsMat, randRow, randCol, randData, failure_info);
}

template <class LayoutType, class ExeSpaceType>
void doAllScalarsCoo2Crs(size_t m, size_t n, int min, int max) {
  doCoo2Crs<float, LayoutType, ExeSpaceType>(m, n, min, max);
  doCoo2Crs<double, LayoutType, ExeSpaceType>(m, n, min, max);
  doCoo2Crs<Kokkos::complex<float>, LayoutType, ExeSpaceType>(m, n, min, max);
  doCoo2Crs<Kokkos::complex<double>, LayoutType, ExeSpaceType>(m, n, min, max);
}

template <class ExeSpaceType>
void doAllLayoutsCoo2Crs(size_t m, size_t n, int min, int max) {
  doAllScalarsCoo2Crs<Kokkos::LayoutLeft, ExeSpaceType>(m, n, min, max);
  doAllScalarsCoo2Crs<Kokkos::LayoutRight, ExeSpaceType>(m, n, min, max);
}

template <class ExeSpaceType>
void doAllCoo2Crs(size_t m, size_t n) {
  int min = 1, max = 10;
  doAllLayoutsCoo2Crs<ExeSpaceType>(m, n, min, max);
}

TEST_F(TestCategory, sparse_coo2crs) {
#if defined(KOKKOS_ENABLE_SYCL)
  if constexpr (std::is_same_v<typename TestDevice::execution_space, Kokkos::Experimental::SYCL>) {
    std::cout << "Not running coo2csr on SYCL execution space" << std::endl;
    return;
  }
#endif

  uint64_t ticks = std::chrono::high_resolution_clock::now().time_since_epoch().count() % UINT32_MAX;
  std::srand(ticks);

  doAllCoo2Crs<TestDevice>(0, 0);

  // Square cases
  for (size_t i = 1; i < 256; i *= 4) {
    size_t dim = (std::rand() % 511) + 1;
    doAllCoo2Crs<TestDevice>(dim, dim);
  }

  // Non-square cases
  for (size_t i = 1; i < 256; i *= 4) {
    size_t m = (std::rand() % 511) + 1;
    size_t n = (std::rand() % 511) + 1;
    while (n == m) n = (std::rand() % 511) + 1;
    doAllCoo2Crs<TestDevice>(m, n);
  }

  RandCooMat<double, Kokkos::LayoutRight, TestDevice> cooMat(2, 2, 2 * 2, 10, 10);
  auto crsMatrix = KokkosSparse::coo2crs(2, 2, cooMat.get_row(), cooMat.get_col(), cooMat.get_data());
  auto cooMatrix = KokkosSparse::crs2coo(crsMatrix);

  check_crs_matrix(crsMatrix, cooMatrix.row(), cooMatrix.col(), cooMatrix.data());
}

TEST_F(TestCategory, sparse_coo2crs_staticMatrix_edgeCases) {
#if defined(KOKKOS_ENABLE_SYCL)
  if constexpr (std::is_same_v<typename TestDevice::execution_space, Kokkos::Experimental::SYCL>) {
    std::cout << "Not running coo2csr on SYCL execution space" << std::endl;
    return;
  }
#endif

  int m = 4;
  int n = 4;
  long long staticRow[16]{0, 1, 3, 2, 3, 2, 2, 2, 0, 0, 0, 1, 2, 0, 3, 0};
  long long staticCol[16]{1, 1, 2, 3, 3, 2, 3, 2, 0, 0, 1, 3, 1, 2, 0, 0};
  float staticData[16]{7.28411, 8.17991, 8.84304, 5.01788, 9.85646, 5.79404, 8.42014, 1.90238,
                       8.24195, 4.39955, 3.2637,  5.4546,  6.51895, 8.09302, 9.36294, 3.44206};
  Kokkos::View<long long *, TestDevice> row("coo row", 16);
  Kokkos::View<long long *, TestDevice> col("coo col", 16);
  Kokkos::View<float *, TestDevice> data("coo data", 16);

  typename Kokkos::View<long long *, TestDevice>::HostMirror row_h = Kokkos::create_mirror_view(row);
  typename Kokkos::View<long long *, TestDevice>::HostMirror col_h = Kokkos::create_mirror_view(col);
  typename Kokkos::View<float *, TestDevice>::HostMirror data_h    = Kokkos::create_mirror_view(data);
  for (int i = 0; i < 16; i++) {
    row_h(i)  = staticRow[i];
    col_h(i)  = staticCol[i];
    data_h(i) = staticData[i];
  }

  Kokkos::deep_copy(row, row_h);
  Kokkos::deep_copy(col, col_h);
  Kokkos::deep_copy(data, data_h);

  // Even partitions with multiple threads
  auto crsMatTs4 = KokkosSparse::coo2crs(m, n, row, col, data);
  check_crs_matrix(crsMatTs4, row_h, col_h, data_h);

  // Even partitions, single thread, fully sparse row
  long long staticRowTs1[16]{0, 3, 0, 2, 2, 3, 0, 3, 2, 0, 0, 0, 0, 3, 3, 0};
  long long staticColTs1[16]{3, 1, 3, 1, 2, 2, 1, 1, 2, 3, 3, 1, 1, 0, 0, 0};
  float staticDataTs1[16]{6.1355,  6.53989, 8.58559, 6.37476, 4.18964, 2.41146, 1.82177, 1.4249,
                          1.52659, 5.50521, 8.0484,  3.98874, 6.74709, 3.35072, 7.81944, 5.83494};
  for (int i = 0; i < 16; i++) {
    row_h(i)  = staticRowTs1[i];
    col_h(i)  = staticColTs1[i];
    data_h(i) = staticDataTs1[i];
  }
  Kokkos::deep_copy(row, row_h);
  Kokkos::deep_copy(col, col_h);
  Kokkos::deep_copy(data, data_h);

  auto crsMatTs1 = KokkosSparse::coo2crs(m, n, row, col, data);
  check_crs_matrix(crsMatTs1, row_h, col_h, data_h);

  // Fully sparse
  for (int i = 0; i < 16; i++) {
    row_h(i) = -staticRowTs1[i];
    col_h(i) = -staticColTs1[i];
  }
  Kokkos::deep_copy(row, row_h);
  Kokkos::deep_copy(col, col_h);

  auto crsMatFsTs1 = KokkosSparse::coo2crs(m, n, row, col, data);
  check_crs_matrix(crsMatFsTs1, row_h, col_h, data);
}
}  // namespace Test
