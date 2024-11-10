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

#ifndef KOKKOSSPARSE_BSR_TO_CRS_IMPL_HPP
#define KOKKOSSPARSE_BSR_TO_CRS_IMPL_HPP

#include "KokkosSparse_BsrMatrix.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

namespace KokkosSparse {

namespace Impl {

/*! \brief Create an equivalent point matrix from a Bsr matrix
    The Crs and Bsr matrix do not have to be on the same device
*/
template <typename Crs, typename Bsr>
Crs bsr_to_crs(const Bsr &bsr) {
  using crs_device_type  = typename Crs::device_type;
  using crs_values_type  = typename Crs::values_type;
  using crs_index_type   = typename Crs::index_type;
  using crs_ordinal_type = typename Crs::non_const_ordinal_type;
  using crs_scalar_type  = typename Crs::non_const_value_type;
  using crs_size_type    = typename Crs::non_const_size_type;

  using crs_row_map_type = Kokkos::View<typename Crs::row_map_type::non_const_data_type, crs_device_type>;
  using bsr_ordinal_type = typename Bsr::non_const_ordinal_type;

  using bsr_size_type = typename Bsr::non_const_size_type;

  // determine what some output matrix parameter will be
  const size_t bs                   = bsr.blockDim();
  const crs_ordinal_type crsNumRows = bsr.numRows() * bs;
  const crs_ordinal_type crsNumCols = bsr.numCols() * bs;
  const crs_size_type crsNnz        = bsr.nnz() * bs * bs;

  // clone Bsr row map to host memory space
  auto bRows = Kokkos::create_mirror_view(bsr.graph.row_map);
  auto bInds = Kokkos::create_mirror_view(bsr.graph.entries);
  auto bVals = Kokkos::create_mirror_view(bsr.values);
  Kokkos::deep_copy(bRows, bsr.graph.row_map);
  Kokkos::deep_copy(bInds, bsr.graph.entries);
  Kokkos::deep_copy(bVals, bsr.values);

  using Entry = std::pair<crs_ordinal_type, crs_scalar_type>;  // {column, value}
  using Row   = std::vector<Entry>;                            // all entries in a row
  std::map<crs_ordinal_type, Row> rows;                        // entries in each row

  // sort entries in a row by column
  auto by_col = [](const Entry &a, const Entry &b) { return a.first < b.first; };

  // Convert BSR data into CRS rows
  for (bsr_ordinal_type bRow = 0; bRow < bsr_ordinal_type(bsr.numRows()); ++bRow) {
    for (bsr_size_type bColIdx = bRows(bRow); bColIdx < bRows(bRow + 1); ++bColIdx) {
      const crs_ordinal_type bCol = bInds(bColIdx);

      // add all points in this block
      for (bsr_size_type lr = 0; lr < bsr_size_type(bs); ++lr) {
        const crs_ordinal_type cRow = bRow * bs + lr;
        for (bsr_size_type lc = 0; lc < bsr_size_type(bs); ++lc) {
          const crs_size_type cvi     = bColIdx * bs * bs + lr * bs + lc;
          const crs_ordinal_type cCol = bCol * bs + lc;
          const crs_scalar_type cVal  = bVals(cvi);
          auto entry                  = std::make_pair(cCol, cVal);

          auto it = rows.find(cRow);
          if (it == rows.end()) {
            Row newRow;
            newRow.push_back(entry);
            rows[cRow] = newRow;
          } else {
            it->second.push_back(entry);
          }
        }
      }
    }
  }

  // device and host views of Crs data
  crs_row_map_type devCrsRows("crs row map", crsNumRows + 1);
  crs_index_type devCrsIdx("crs columns", crsNnz);
  crs_values_type devCrsVals("crs values", crsNnz);
  auto hostCrsRows = Kokkos::create_mirror_view(devCrsRows);
  auto hostCrsIdx  = Kokkos::create_mirror_view(devCrsIdx);
  auto hostCrsVals = Kokkos::create_mirror_view(devCrsVals);

  // convert to Crs format
  crs_ordinal_type iRowMap = 0;
  crs_size_type nentries   = 0;
  for (auto &kv : rows) {                     // iterating through rows in order
    const crs_ordinal_type &row = kv.first;   // block's position
    Row &entries                = kv.second;  // non-zeros in the block

    // update row map if we've moved to a new row
    for (; iRowMap < row; ++iRowMap) {
      hostCrsRows(iRowMap + 1) = nentries;  // row ends at entries so far
    }

    // make sure crs points in each row are sorted by column
    std::sort(entries.begin(), entries.end(), by_col);

    // add columns and values to Crs data
    for (size_t i = 0; i < entries.size(); ++i, ++nentries) {
      hostCrsIdx(nentries)  = entries[i].first;
      hostCrsVals(nentries) = entries[i].second;
    }
  }
  // complete row map if last blocks are empty
  for (; iRowMap < crsNumRows; ++iRowMap) {
    hostCrsRows(iRowMap + 1) = nentries;
  }

  // move to device
  Kokkos::deep_copy(devCrsRows, hostCrsRows);
  Kokkos::deep_copy(devCrsIdx, hostCrsIdx);
  Kokkos::deep_copy(devCrsVals, hostCrsVals);

  // construct the resulting Crs matrix
  Crs crs("", crsNumRows, crsNumCols, crsNnz, devCrsVals, devCrsRows, devCrsIdx);
  return crs;
}  // bsr_to_crs

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_BSR_TO_CRS_IMPL_HPP