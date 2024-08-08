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
#ifndef KOKKOSSPARSE_CRS_DETECT_BLOCK_SIZE_HPP
#define KOKKOSSPARSE_CRS_DETECT_BLOCK_SIZE_HPP

#include <map>

#include <Kokkos_Core.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_Utils.hpp"

/*! \file KokkosSparse_crs_detect_block_size.hpp

    \brief A utility function for detecting the block size in a CrsMatrix. Not
   for performance-sensitive use.
*/

namespace KokkosSparse::Impl {

/**
 * \class BlockPopulations
 * \brief A class to store population counts of blocks in a CrsMatrix
 */
class BlockPopulations {
 public:
  /**
   * \brief Constructor for BlockPopulations
   * \param sz The block size
   */
  BlockPopulations(size_t sz) : sz_(sz) {}

  /**
   * \brief Add a point to the corresponding block
   * \param r The row index of the point
   * \param c The column index of the point
   */
  void add(size_t r, size_t c) {
    auto key = std::make_pair(r / sz_, c / sz_);
    auto it  = blocks_.find(key);
    if (it == blocks_.end()) {
      blocks_.insert(std::make_pair(key, 1));
    } else {
      ++(it->second);
    }
  }

  /**
   * \brief Check if all blocks are dense
   * \return True if all blocks have a count equal to the block size squared
   */
  bool all_dense() const {
    for (const auto &kv : blocks_) {
      if (kv.second < sz_ * sz_) {
        return false;
      }
    }
    return true;
  }

 private:
  std::map<std::pair<size_t, size_t>, size_t> blocks_; /**< A map of block coordinates to their population counts */
  size_t sz_;                                          /**< The block size */
};

/**
 * @brief Detects the largest block size that yields only dense blocks in a
 CrsMatrix
 *
 * @tparam Crs The type of the CRS matrix.
 * @param crs The CRS matrix to detect the block size for.
 * @return The largest block size that results in completely dense blocks
    The smallest valid block size is 1
    Since blocks must be dense, sqrt(nnz), num rows, num cols, and min nnz/row
 among non-empty rows are all easy upper bounds of the block size.
 Block sizes are tested from 1 to the minimum of the above.
 The matrix dimensions must divide  evenly into a trial block size (otherwise a
 block would not be full). Furthermore, if a block size of N is not dense, any
 multiple of N will also not be dense, and can be skipped. This is because
 blocks of 2N contain blocks of N, at least one of which is already known not to
 be dense. In practice, this ends up testing only small composite factors and
 all prime factors up to the upper bound.
*/
template <typename Crs>
size_t detect_block_size(const Crs &crs) {
  using ordinal_type = typename Crs::ordinal_type;

  // copy matrix data to host
  auto rs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), crs.graph.row_map);
  auto cs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), crs.graph.entries);

  // upper bound is minimum of sqrt(nnz), numRows, numCols,
  // and smallest non-empty row
  size_t upperBound = std::sqrt(double(crs.nnz()));
  upperBound        = std::min(upperBound, size_t(crs.numRows()));
  upperBound        = std::min(upperBound, size_t(crs.numCols()));
  for (size_t i = 1; i < rs.size(); ++i) {
    size_t rowLen = rs(i) - rs(i - 1);
    if (rowLen > 0) {
      upperBound = std::min(upperBound, rowLen);
    }
  }

  // trial blocks sizes that didn't work out
  std::vector<size_t> rejectedSizes;

  size_t largestBlockSize = 1;  // always a valid block size
  for (size_t trialSize = 2; trialSize <= upperBound; ++trialSize) {
    // trial size must be factor of rows / cols
    if ((crs.numRows() % trialSize) || (crs.numCols() % trialSize)) {
      rejectedSizes.push_back(trialSize);
      continue;
    }

    // trial size must not be a multiple of previously-rejected size
    if (std::any_of(rejectedSizes.begin(), rejectedSizes.end(), [&](size_t f) { return trialSize % f == 0; })) {
      rejectedSizes.push_back(trialSize);
      continue;
    }

    // count the population of all blocks
    BlockPopulations pops(trialSize);
    for (ordinal_type row = 0; row < crs.numRows(); ++row) {
      for (size_t ci = rs(row); ci < rs(row + 1); ++ci) {
        ordinal_type col = cs(ci);
        pops.add(row, col);
      }
    }

    // if all blocks are dense, this is the largest one so far
    if (pops.all_dense()) {
      largestBlockSize = trialSize;
    } else {
      rejectedSizes.push_back(trialSize);
    }
  }
  return largestBlockSize;
}

}  // namespace KokkosSparse::Impl

#endif  // KOKKOSSPARSE_CRS_DETECT_BLOCK_SIZE_HPP
