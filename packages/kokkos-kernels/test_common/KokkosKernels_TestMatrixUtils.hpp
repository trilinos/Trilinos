// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_TESTMATRIXUTILS_HPP
#define KOKKOSKERNELS_TESTMATRIXUTILS_HPP

#include <random>
#include <vector>

#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_ArithTraits.hpp"
#include "KokkosBatched_Vector.hpp"
// Make this include-able from all subdirectories

namespace Test {

template <typename scalar_t, typename lno_t, typename size_type, typename device, typename crsMat_t>
crsMat_t symmetrize(crsMat_t A) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  auto host_rowmap  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
  auto host_entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);
  auto host_values  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
  lno_t numRows     = A.numRows();
  // symmetrize as input_mat + input_mat^T, to still have a diagonally dominant
  // matrix
  typedef std::map<lno_t, scalar_t> Row;
  std::vector<Row> symRows(numRows);
  for (lno_t r = 0; r < numRows; r++) {
    auto& row = symRows[r];
    for (size_type i = host_rowmap(r); i < host_rowmap(r + 1); i++) {
      lno_t c   = host_entries(i);
      auto& col = symRows[c];
      auto it   = row.find(c);
      if (it == row.end())
        row[c] = host_values(i);
      else
        row[c] += host_values(i);
      it = col.find(r);
      if (it == col.end())
        col[r] = host_values(i);
      else
        col[r] += host_values(i);
    }
  }
  // Count entries
  Kokkos::View<size_type*, Kokkos::LayoutLeft, Kokkos::HostSpace> new_host_rowmap("Rowmap", numRows + 1);
  size_t accum = 0;
  for (lno_t r = 0; r <= numRows; r++) {
    new_host_rowmap(r) = accum;
    if (r < numRows) accum += symRows[r].size();
  }
  // Allocate new entries/values
  Kokkos::View<lno_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> new_host_entries("Entries", accum);
  Kokkos::View<scalar_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> new_host_values("Values", accum);
  for (lno_t r = 0; r < numRows; r++) {
    auto rowIt = symRows[r].begin();
    for (size_type i = new_host_rowmap(r); i < new_host_rowmap(r + 1); i++) {
      new_host_entries(i) = rowIt->first;
      new_host_values(i)  = rowIt->second;
      rowIt++;
    }
  }
  lno_view_t new_rowmap("Rowmap", numRows + 1);
  lno_nnz_view_t new_entries("Entries", accum);
  scalar_view_t new_values("Values", accum);
  Kokkos::deep_copy(new_rowmap, new_host_rowmap);
  Kokkos::deep_copy(new_entries, new_host_entries);
  Kokkos::deep_copy(new_values, new_host_values);
  return crsMat_t("SymA", numRows, numRows, accum, new_values, new_rowmap, new_entries);
}

/**
 * Utilities for defining small matrix fixtures as 2D std::vectors, and
 * converting them to and from compressed sparse formats (CRS/CSC) for use in
 * unit tests.
 *
 * Fixtures are defined as row-major 2D vectors of scalars, where zeros
 * represent structural zeros (absent entries). Use KEEP_ZERO() as a sentinel
 * value when an explicit zero entry must be stored in the sparse matrix.
 *
 * Example — define a 4x4 CRS matrix and verify its contents:
 *
 *   using scalar_t  = double;
 *   using size_type = int;
 *   using rowmap_t  = Kokkos::View<size_type*>;
 *   using entries_t = Kokkos::View<size_type*>;
 *   using values_t  = Kokkos::View<scalar_t*>;
 *
 *   std::vector<std::vector<scalar_t>> fixture = {
 *     { 1.0, 2.0, 0.0, 0.0 },
 *     { 0.0, 3.0, 4.0, 0.0 },
 *     { 0.0, 0.0, 5.0, 6.0 },
 *     { 7.0, 0.0, 0.0, 8.0 },
 *   };
 *
 *   rowmap_t  row_map("row_map", 0);
 *   entries_t entries("entries", 0);
 *   values_t  values("values",  0);
 *   compress_matrix(row_map, entries, values, fixture);
 *
 *   // ... apply kernel under test ...
 *
 *   check_matrix("result", row_map, entries, values, fixture);
 *
 * To build a CSC matrix instead, pass true as the CSC template argument:
 *
 *   compress_matrix<true>(col_map, entries, values, fixture);
 */

template <typename scalar_t>
scalar_t KEEP_ZERO() {
  return scalar_t(-9999.0);
}

template <bool CSC = false, typename MapT, typename EntriesT, typename ValuesT>
void compress_matrix(MapT& map, EntriesT& entries, ValuesT& values,
                     const std::vector<std::vector<typename ValuesT::non_const_value_type>>& fixture) {
  using size_type = typename MapT::non_const_value_type;
  using scalar_t  = typename ValuesT::non_const_value_type;

  const scalar_t ZERO = scalar_t(0);

  const size_type nrows = fixture.size();
  const size_type ncols = fixture[0].size();

  // Count fixture nnz's
  size_type nnz = 0;
  for (size_type row_idx = 0; row_idx < nrows; ++row_idx) {
    for (size_type col_idx = 0; col_idx < nrows; ++col_idx) {
      if (fixture[row_idx][col_idx] != ZERO) {
        ++nnz;
      }
    }
  }

  // Allocate device CRS views
  Kokkos::resize(map, (CSC ? ncols : nrows) + 1);
  Kokkos::resize(entries, nnz);
  Kokkos::resize(values, nnz);

  // Create host mirror views for CRS
  auto hmap     = Kokkos::create_mirror_view(map);
  auto hentries = Kokkos::create_mirror_view(entries);
  auto hvalues  = Kokkos::create_mirror_view(values);

  // Compress into CRS (host views)
  size_type curr_nnz = 0;

  const size_type num_outer = (CSC ? ncols : nrows);
  const size_type num_inner = (CSC ? nrows : ncols);
  for (size_type outer_idx = 0; outer_idx < num_outer; ++outer_idx) {
    for (size_type inner_idx = 0; inner_idx < num_inner; ++inner_idx) {
      const size_type row = CSC ? inner_idx : outer_idx;
      const size_type col = CSC ? outer_idx : inner_idx;
      const auto val      = fixture[row][col];
      if (val != ZERO) {
        hentries(curr_nnz) = inner_idx;
        hvalues(curr_nnz)  = val == KEEP_ZERO<scalar_t>() ? ZERO : val;
        ++curr_nnz;
      }
      hmap(outer_idx + 1) = curr_nnz;
    }
  }

  // Copy host CRS views to device CRS views
  Kokkos::deep_copy(map, hmap);
  Kokkos::deep_copy(entries, hentries);
  Kokkos::deep_copy(values, hvalues);
}

template <bool CSC = false, typename RowMapT, typename EntriesT, typename ValuesT>
std::vector<std::vector<typename ValuesT::non_const_value_type>> decompress_matrix(const RowMapT& row_map,
                                                                                   const EntriesT& entries,
                                                                                   const ValuesT& values) {
  using size_type = typename RowMapT::non_const_value_type;
  using scalar_t  = typename ValuesT::non_const_value_type;

  const scalar_t ZERO = scalar_t(0);

  const size_type nrows = row_map.size() - 1;
  std::vector<std::vector<scalar_t>> result;
  result.resize(nrows);
  for (auto& row : result) {
    row.resize(nrows, ZERO);
  }

  auto hrow_map = Kokkos::create_mirror_view(row_map);
  auto hentries = Kokkos::create_mirror_view(entries);
  auto hvalues  = Kokkos::create_mirror_view(values);
  Kokkos::deep_copy(hrow_map, row_map);
  Kokkos::deep_copy(hentries, entries);
  Kokkos::deep_copy(hvalues, values);

  for (size_type row_idx = 0; row_idx < nrows; ++row_idx) {
    const size_type row_nnz_begin = hrow_map(row_idx);
    const size_type row_nnz_end   = hrow_map(row_idx + 1);
    for (size_type row_nnz = row_nnz_begin; row_nnz < row_nnz_end; ++row_nnz) {
      const auto col_idx   = hentries(row_nnz);
      const scalar_t value = hvalues(row_nnz);
      if (CSC) {
        result[col_idx][row_idx] = value;
      } else {
        result[row_idx][col_idx] = value;
      }
    }
  }

  return result;
}

template <typename RowMapT, typename EntriesT, typename ValuesT>
std::vector<std::vector<typename ValuesT::non_const_value_type>> decompress_matrix(
    const RowMapT& row_map, const EntriesT& entries, const ValuesT& values,
    typename RowMapT::const_value_type block_size) {
  using size_type = typename RowMapT::non_const_value_type;
  using scalar_t  = typename ValuesT::non_const_value_type;

  const scalar_t ZERO = scalar_t(0);

  const size_type nbrows      = row_map.extent(0) - 1;
  const size_type nrows       = nbrows * block_size;
  const size_type block_items = block_size * block_size;
  std::vector<std::vector<scalar_t>> result;
  result.resize(nrows);
  for (auto& row : result) {
    row.resize(nrows, ZERO);
  }

  auto hrow_map = Kokkos::create_mirror_view(row_map);
  auto hentries = Kokkos::create_mirror_view(entries);
  auto hvalues  = Kokkos::create_mirror_view(values);
  Kokkos::deep_copy(hrow_map, row_map);
  Kokkos::deep_copy(hentries, entries);
  Kokkos::deep_copy(hvalues, values);

  for (size_type row_idx = 0; row_idx < nbrows; ++row_idx) {
    const size_type row_nnz_begin = hrow_map(row_idx);
    const size_type row_nnz_end   = hrow_map(row_idx + 1);
    for (size_type row_nnz = row_nnz_begin; row_nnz < row_nnz_end; ++row_nnz) {
      const auto col_idx = hentries(row_nnz);
      for (size_type i = 0; i < block_size; ++i) {
        const size_type unc_row_idx = row_idx * block_size + i;
        for (size_type j = 0; j < block_size; ++j) {
          const size_type unc_col_idx      = col_idx * block_size + j;
          result[unc_row_idx][unc_col_idx] = hvalues(row_nnz * block_items + i * block_size + j);
        }
      }
    }
  }

  return result;
}

template <typename RowMapT, typename EntriesT, typename ValuesT>
void check_matrix(const std::string& name, const RowMapT& row_map, const EntriesT& entries, const ValuesT& values,
                  const std::vector<std::vector<typename ValuesT::non_const_value_type>>& expected) {
  using size_type = typename RowMapT::non_const_value_type;

  const auto decompressed_mtx = decompress_matrix(row_map, entries, values);

  const size_type nrows = row_map.size() - 1;
  for (size_type row_idx = 0; row_idx < nrows; ++row_idx) {
    for (size_type col_idx = 0; col_idx < nrows; ++col_idx) {
      EXPECT_NEAR(expected[row_idx][col_idx], decompressed_mtx[row_idx][col_idx], 0.01)
          << "Failed check is: " << name << "[" << row_idx << "][" << col_idx << "]";
    }
  }
}

template <typename scalar_t>
void print_matrix(const std::vector<std::vector<scalar_t>>& matrix) {
  for (const auto& row : matrix) {
    for (const auto& item : row) {
      std::printf("%.5f ", item);
    }
    std::cout << std::endl;
  }
}

/// Fill a 2D view from a row-major 2D vector fixture.
/// Works for both host and device views: internally creates a host mirror,
/// fills it, then deep-copies to the view.
template <typename ViewType>
void fill_view_from_fixture(ViewType& view,
                            const std::vector<std::vector<typename ViewType::non_const_value_type>>& fixture) {
  auto host_view = Kokkos::create_mirror_view(view);
  for (size_t i = 0; i < fixture.size(); ++i)
    for (size_t j = 0; j < fixture[i].size(); ++j) host_view(i, j) = fixture[i][j];
  Kokkos::deep_copy(view, host_view);
}

}  // namespace Test
#endif
