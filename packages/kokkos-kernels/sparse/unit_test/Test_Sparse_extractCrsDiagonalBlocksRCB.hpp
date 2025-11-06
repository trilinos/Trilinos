// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "KokkosSparse_Utils.hpp"
#include "KokkosKernels_TestUtils.hpp"

namespace Test {
template <typename coors_view_t>
typename coors_view_t::value_type generate_3d_coordinates_for_sparse_rows(int n_pts_per_dim,
                                                                          coors_view_t &coordinates) {
  using scalar_t = typename coors_view_t::value_type;

  scalar_t x_max = 1.1;
  scalar_t x_min = 0;
  scalar_t y_max = 1.05;
  scalar_t y_min = 0;
  scalar_t z_max = 1.0;
  scalar_t z_min = 0;

  int n_spaces = n_pts_per_dim - 1;

  scalar_t dx = (x_max - x_min) / n_spaces;
  scalar_t dy = (y_max - y_min) / n_spaces;
  scalar_t dz = (z_max - z_min) / n_spaces;

  auto h_coordinates = Kokkos::create_mirror_view(coordinates);

  int cnt = 0;
  for (int ii = 0; ii < n_pts_per_dim; ii++) {  // z
    scalar_t z = z_min + ii * dz;
    for (int jj = 0; jj < n_pts_per_dim; jj++) {  // y
      scalar_t y = y_min + jj * dy;
      for (int kk = 0; kk < n_pts_per_dim; kk++) {  // x
        scalar_t x            = x_min + kk * dx;
        h_coordinates(cnt, 0) = x;
        h_coordinates(cnt, 1) = y;
        h_coordinates(cnt, 2) = z;
        cnt++;
      }
    }
  }

  Kokkos::deep_copy(coordinates, h_coordinates);

  return dx;
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void run_test_extract_diagonal_blocks_rcb(lno_t n_pts_per_dim, lno_t nblocks) {
  using RowMapType       = Kokkos::View<size_type *, device>;
  using EntriesType      = Kokkos::View<lno_t *, device>;
  using ValuesType       = Kokkos::View<scalar_t *, device>;
  using magnitude_t      = typename KokkosKernels::ArithTraits<scalar_t>::mag_type;
  using CoorsViewType    = Kokkos::View<magnitude_t **, device>;
  using PermViewType     = Kokkos::View<lno_t *, device>;
  using CoorsViewType_hm = typename CoorsViewType::host_mirror_type;
  using PermViewType_hm  = typename PermViewType::host_mirror_type;
  using crsMat_t         = KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type>;

  crsMat_t A;
  std::vector<crsMat_t> DiagBlks(nblocks);

  // Generate coordinates
  lno_t n_coordinates = n_pts_per_dim * n_pts_per_dim * n_pts_per_dim;
  CoorsViewType coordinates(Kokkos::view_alloc(Kokkos::WithoutInitializing, "coordinates"), n_coordinates, 3);
  CoorsViewType_hm h_coordinates = Kokkos::create_mirror(coordinates);
  magnitude_t dx = generate_3d_coordinates_for_sparse_rows<CoorsViewType_hm>(n_pts_per_dim, h_coordinates);
  Kokkos::deep_copy(coordinates, h_coordinates);

  // Generate test matrix consisting of near interactions calculated based on coordinates
  lno_t nrows = n_coordinates;
  RowMapType row_map(Kokkos::view_alloc(Kokkos::WithoutInitializing, "row_map"), nrows + 1);
  auto h_row_map = Kokkos::create_mirror_view(row_map);

  std::vector<std::pair<lno_t, lno_t>> near_indices;
  lno_t nnz = 0;
  for (lno_t i = 0; i < nrows; i++) {
    h_row_map(i) = nnz;
    for (lno_t j = 0; j < nrows; j++) {
      magnitude_t distance = sqrt(pow(h_coordinates(i, 0) - h_coordinates(j, 0), 2.0) +
                                  pow(h_coordinates(i, 1) - h_coordinates(j, 1), 2.0) +
                                  pow(h_coordinates(i, 2) - h_coordinates(j, 2), 2.0));
      if (distance <= dx) {
        nnz++;
        near_indices.push_back(std::make_pair(i, j));
      }
    }
  }
  h_row_map(nrows) = nnz;

  EntriesType entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entries"), nnz);
  ValuesType values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "values"), nnz);
  auto h_entries = Kokkos::create_mirror_view(entries);
  auto h_values  = Kokkos::create_mirror_view(values);

  for (lno_t ii = 0; ii < nnz; ii++) {
    lno_t i       = near_indices[ii].first;
    lno_t j       = near_indices[ii].second;
    h_entries(ii) = j;
    h_values(ii)  = i * 10 + j + 1;
  }

  // Copy from host to device
  Kokkos::deep_copy(row_map, h_row_map);
  Kokkos::deep_copy(entries, h_entries);
  Kokkos::deep_copy(values, h_values);

  // Construct a CRS matrix
  A = crsMat_t("CrsMatrix", nrows, nrows, nnz, values, row_map, entries);

  // Extract diagonal blocks
  PermViewType perm_rcb(Kokkos::view_alloc(Kokkos::WithoutInitializing, "perm_rcb"), n_coordinates);
  KokkosSparse::Impl::kk_extract_diagonal_blocks_crsmatrix_with_rcb_sequential(A, coordinates, DiagBlks, perm_rcb);

  // Checking results
  lno_t numRows = 0;
  lno_t numCols = 0;
  for (lno_t i = 0; i < nblocks; i++) {
    numRows += DiagBlks[i].numRows();
    numCols += DiagBlks[i].numCols();
  }

  ASSERT_EQ(numRows, nrows);
  ASSERT_EQ(numCols, nrows);

  lno_t n_levels = static_cast<lno_t>(std::log2(static_cast<double>(nblocks)) + 1);
  PermViewType_hm perm_rcb_ref(Kokkos::view_alloc(Kokkos::WithoutInitializing, "perm_rcb_ref"), n_coordinates);
  PermViewType_hm reverse_perm_rcb_ref(Kokkos::view_alloc(Kokkos::WithoutInitializing, "reverse_perm_rcb_ref"),
                                       n_coordinates);
  std::vector<lno_t> partition_sizes = KokkosGraph::Experimental::recursive_coordinate_bisection(
      h_coordinates, perm_rcb_ref, reverse_perm_rcb_ref, n_levels);

  std::map<lno_t, scalar_t> colIdx_Value_rcb;

  lno_t blk_size;
  lno_t blk_start = 0;
  for (lno_t i = 0; i < nblocks; i++) {  // block loop
    auto h_row_map_diagblk = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DiagBlks[i].graph.row_map);
    auto h_entries_diagblk = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DiagBlks[i].graph.entries);
    auto h_values_diagblk  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), DiagBlks[i].values);

    ASSERT_EQ(static_cast<lno_t>(DiagBlks[i].numRows()), partition_sizes[i]);

    blk_size          = partition_sizes[i];
    size_type blk_nnz = 0;
    for (lno_t ii = 0; ii < blk_size; ii++) {  // row loop in each block
      ASSERT_EQ(h_row_map_diagblk(ii), blk_nnz);
      colIdx_Value_rcb.clear();
      lno_t origRow = reverse_perm_rcb_ref(blk_start + ii);  // get the original row idx of the reordered row idx, ii
      for (size_type j = h_row_map(origRow); j < h_row_map(origRow + 1); j++) {
        lno_t origEi   = h_entries(j);
        scalar_t origV = h_values(j);
        lno_t Ei       = perm_rcb_ref(origEi);  // get the reordered col idx of the
                                                // original col idx, origEi
        colIdx_Value_rcb[Ei] = origV;
      }
      for (typename std::map<lno_t, scalar_t>::iterator it = colIdx_Value_rcb.begin(); it != colIdx_Value_rcb.end();
           ++it) {
        if ((it->first >= blk_start) && (it->first < (blk_start + blk_size))) {
          ASSERT_EQ(h_entries_diagblk(blk_nnz) + blk_start, it->first);
          ASSERT_EQ(h_values_diagblk(blk_nnz), it->second);
          blk_nnz++;
        }
      }
    }  // row loop in each block
    blk_start += DiagBlks[i].numCols();
  }  // block loop
}
}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_extract_diagonal_blocks_rcb() {
  Test::run_test_extract_diagonal_blocks_rcb<scalar_t, lno_t, size_type, device>(5, 2);
  Test::run_test_extract_diagonal_blocks_rcb<scalar_t, lno_t, size_type, device>(5, 4);
  Test::run_test_extract_diagonal_blocks_rcb<scalar_t, lno_t, size_type, device>(5, 8);
  Test::run_test_extract_diagonal_blocks_rcb<scalar_t, lno_t, size_type, device>(5, 16);
  Test::run_test_extract_diagonal_blocks_rcb<scalar_t, lno_t, size_type, device>(9, 2);
  Test::run_test_extract_diagonal_blocks_rcb<scalar_t, lno_t, size_type, device>(9, 4);
  Test::run_test_extract_diagonal_blocks_rcb<scalar_t, lno_t, size_type, device>(9, 8);
  Test::run_test_extract_diagonal_blocks_rcb<scalar_t, lno_t, size_type, device>(9, 16);
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                           \
  TEST_F(TestCategory, sparse##_##extract_diagonal_blocks_rcb##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_extract_diagonal_blocks_rcb<SCALAR, ORDINAL, OFFSET, DEVICE>();                                      \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
