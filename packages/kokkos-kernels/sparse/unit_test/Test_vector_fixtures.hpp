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

#ifndef _TEST_VECTOR_FIXTURES_HPP
#define _TEST_VECTOR_FIXTURES_HPP

#include <Kokkos_Core.hpp>

#include <vector>

/**
 * API for working with 2D vectors of small matrices for testing.
 */

namespace Test {

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

}  // namespace Test

#endif  // _TEST_VECTOR_FIXTURES_HPP
