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

#include <sstream>
#include <iostream>
#include <iomanip>

#include "Kokkos_Core.hpp"

#include "KokkosKernels_default_types.hpp"
#include "KokkosSparse_BsrMatrix.hpp"

using Scalar  = default_scalar;
using Ordinal = default_lno_t;
using Offset  = default_size_type;
using Layout  = default_layout;

template <class bsrmatrix_type>
struct bsr_fill {
  bsrmatrix_type bsr_mat;

  bsr_fill(bsrmatrix_type bsr_mat_) : bsr_mat(bsr_mat_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& rowIdx) const {
    if (rowIdx == 0) {  // Left boundary condition
      auto block_tmp  = bsr_mat.unmanaged_block(0);
      block_tmp(0, 0) = 1.0;
      block_tmp(0, 1) = 0.0;
      block_tmp(1, 0) = 0.0;
      block_tmp(1, 1) = 1.0;
    } else if (rowIdx == bsr_mat.numRows() - 1) {  // Right boundary condition
      auto block_tmp  = bsr_mat.unmanaged_block(bsr_mat.graph.row_map(rowIdx) + 1);
      block_tmp(0, 0) = 1.0;
      block_tmp(1, 1) = 1.0;
    } else {
      auto block_tmp  = bsr_mat.unmanaged_block(bsr_mat.graph.row_map(rowIdx));
      block_tmp(0, 0) = -1.0;
      block_tmp(0, 1) = -1.0 / 2.0;
      block_tmp(1, 0) = 0.0;
      block_tmp(1, 1) = -1.0;

      block_tmp       = bsr_mat.unmanaged_block(bsr_mat.graph.row_map(rowIdx) + 1);
      block_tmp(0, 0) = 2.0;
      block_tmp(0, 1) = 0.0;
      block_tmp(1, 0) = 0.0;
      block_tmp(1, 1) = 2.0;

      block_tmp       = bsr_mat.unmanaged_block(bsr_mat.graph.row_map(rowIdx) + 2);
      block_tmp(0, 0) = -1.0;
      block_tmp(0, 1) = 1.0 / 2.0;
      block_tmp(1, 0) = 0.0;
      block_tmp(1, 1) = -1.0;
    }
  }
};

template <class bsrmatrix_type, class diag_blocks_type>
struct diagonal_extractor {
  using graph_type     = typename bsrmatrix_type::staticcrsgraph_type;
  using row_map_type   = typename graph_type::row_map_type;
  using entries_type   = typename graph_type::entries_type;
  using bsr_block_type = typename bsrmatrix_type::block_type;

  bsrmatrix_type bsr_mat;
  row_map_type row_map;
  entries_type entries;
  diag_blocks_type diag_blocks;

  diagonal_extractor(bsrmatrix_type bsr_mat_, diag_blocks_type diag_blocks_)
      : bsr_mat(bsr_mat_),
        row_map(bsr_mat_.graph.row_map),
        entries(bsr_mat_.graph.entries),
        diag_blocks(diag_blocks_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& rowIdx) const {
    for (Offset entryIdx = row_map(rowIdx); entryIdx < row_map(rowIdx + 1); ++entryIdx) {
      if (entries(entryIdx) == rowIdx) {
        bsr_block_type bsr_diag_block = bsr_mat.unmanaged_block(entryIdx);
        for (int i = 0; i < bsr_mat.blockDim(); ++i) {
          for (int j = 0; j < bsr_mat.blockDim(); ++j) {
            diag_blocks(rowIdx, i, j) = bsr_diag_block(i, j);
          }
        }
      }
    }
  }
};

int main(int argc, char* argv[]) {
  using device_type =
      typename Kokkos::Device<Kokkos::DefaultExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory_space>;
  using bsrmatrix_type = typename KokkosSparse::Experimental::BsrMatrix<Scalar, Ordinal, device_type, void, Offset>;
  using graph_type     = typename bsrmatrix_type::staticcrsgraph_type;
  using row_map_type   = typename graph_type::row_map_type;
  using entries_type   = typename graph_type::entries_type;

  Kokkos::initialize(argc, argv);
  {
    //
    // We will create a 1D discretization for the coupled thermo-elastic
    // diffusion
    //
    //    -\div(EA \grad_s(u) - \alpha(T-T0)I) = f_u
    //                        -\kappa\Delta(T) = f_T
    //
    // The problem is discretized using finite differences as follows:
    //    \frac{d^2 u}{dx^2}\approx \frac{u_{i+1}-2u_i+u_{i-1}}{h_x^2}
    //    \frac{dT}{dx}\approx\frac{T_{i+1}-T_{i-1}}{2h_x}
    //    \frac{d^2T}{dx^2}\approx\frac{T_{i+1}-2T_i+T_{i-1}}{h_x^2}
    //
    // This leads to the combined stencil (assuming all unit coefficients):
    //
    // [-1  1/2] [2 0] [-1  -1/2]
    // [ 0   -1] [0 2] [ 0    -1]
    //
    // First the graph for the mesh will be constructed.
    // Second a BsrMatrix will be constructed from the graph
    // Third the values of the BsrMatrix will be filled.

    constexpr Ordinal blockSize = 2;
    constexpr Ordinal numRows   = 10;
    constexpr Offset numNNZ     = 3 * numRows - 2;
    bsrmatrix_type bsr_mat;

    {
      typename row_map_type::non_const_type row_map(Kokkos::view_alloc(Kokkos::WithoutInitializing, "row pointers"),
                                                    numRows + 1);
      typename entries_type::non_const_type entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "column indices"),
                                                    numNNZ);
      typename row_map_type::HostMirror row_map_h = Kokkos::create_mirror_view(row_map);
      typename entries_type::HostMirror entries_h = Kokkos::create_mirror_view(entries);

      // First Step: build the CrsGraph
      {
        // Build the row pointers and store numNNZ

        row_map_h(0) = 0;
        for (Ordinal rowIdx = 0; rowIdx < numRows; ++rowIdx) {
          if (rowIdx == 0) {
            row_map_h(rowIdx + 1) = row_map_h(rowIdx) + 2;

            entries_h(row_map_h(rowIdx))     = rowIdx;
            entries_h(row_map_h(rowIdx) + 1) = rowIdx + 1;
          } else if (rowIdx == numRows - 1) {
            row_map_h(rowIdx + 1) = row_map_h(rowIdx) + 2;

            entries_h(row_map_h(rowIdx))     = rowIdx - 1;
            entries_h(row_map_h(rowIdx) + 1) = rowIdx;
          } else {
            row_map_h(rowIdx + 1) = row_map_h(rowIdx) + 3;

            entries_h(row_map_h(rowIdx))     = rowIdx - 1;
            entries_h(row_map_h(rowIdx) + 1) = rowIdx;
            entries_h(row_map_h(rowIdx) + 2) = rowIdx + 1;
          }
        }

        if (row_map_h(numRows) != numNNZ) {
          std::ostringstream error_msg;
          error_msg << "error: row_map(numRows) != numNNZ, row_map_h(numRows)=" << row_map_h(numRows)
                    << ", numNNZ=" << numNNZ;
          throw std::runtime_error(error_msg.str());
        }
        Kokkos::deep_copy(row_map, row_map_h);
        Kokkos::deep_copy(entries, entries_h);
      }

      graph_type myGraph(entries, row_map);

      // Second Step: build the BsrMatrix from graph and block size
      bsr_mat = bsrmatrix_type("block matrix", myGraph, blockSize);

      bsr_fill fillFunctor(bsr_mat);
      Kokkos::parallel_for(Kokkos::RangePolicy<int>(0, numRows), fillFunctor);

      std::cout << "BsrMatrix graph: " << std::endl;
      for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
        std::cout << "  [";
        for (int colIdx = 0; colIdx < entries_h(row_map_h(rowIdx)); ++colIdx) {
          std::cout << " ";
        }
        std::cout << "*";
        for (Offset entryIdx = row_map_h(rowIdx); entryIdx < row_map_h(rowIdx + 1) - 1; ++entryIdx) {
          for (int colIdx = entries_h(entryIdx) + 1; colIdx < entries_h(entryIdx + 1); ++colIdx) {
            std::cout << " ";
          }
          std::cout << "*";
        }
        for (int colIdx = entries_h(row_map_h(rowIdx + 1) - 1) + 1; colIdx < numRows; ++colIdx) {
          std::cout << " ";
        }
        std::cout << "]" << std::endl;
      }
    }

    // Extract diagonal block and store them in a rank-3 view
    using diag_blocks_type = Kokkos::View<Scalar***, typename bsrmatrix_type::block_layout_type, device_type>;
    diag_blocks_type diag_blocks("diagonal blocks", numRows, blockSize, blockSize);
    diagonal_extractor myFunc(bsr_mat, diag_blocks);
    Kokkos::parallel_for(Kokkos::RangePolicy<int>(0, numRows), myFunc);

    auto diag_blocks_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, diag_blocks);

    std::cout << "\nBsrMatrix diagonal blocks: " << std::endl;
    for (int blockId = 0; blockId < diag_blocks_h.extent_int(0); ++blockId) {
      std::cout << "  [" << diag_blocks_h(blockId, 0, 0) << ", " << diag_blocks_h(blockId, 0, 1) << "]" << std::endl;
      std::cout << "  [" << diag_blocks_h(blockId, 1, 0) << ", " << diag_blocks_h(blockId, 1, 1) << "]\n" << std::endl;
    }
  }
  Kokkos::finalize();

  return 0;
}
