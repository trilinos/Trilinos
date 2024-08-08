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

#ifndef KOKKOSKERNELS_TEST_STRUCTURE_MATRIX_HPP
#define KOKKOSKERNELS_TEST_STRUCTURE_MATRIX_HPP

#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_Error.hpp"
namespace Test {

enum { FD, FE };

template <class CrsMatrix_t>
struct fill_1D_matrix_functor {
  // Define types used by the CrsMatrix
  typedef typename CrsMatrix_t::execution_space execution_space;
  typedef typename CrsMatrix_t::row_map_type::non_const_type row_map_view_t;
  typedef typename CrsMatrix_t::index_type::non_const_type cols_view_t;
  typedef typename CrsMatrix_t::values_type::non_const_type scalar_view_t;
  typedef typename CrsMatrix_t::non_const_size_type size_type;
  typedef typename CrsMatrix_t::non_const_ordinal_type ordinal_type;

  // Dispatch tags
  struct interiorTag {};
  struct exteriorTag {};

  // Internal variables and temporaries
  const ordinal_type numNodes;
  const int leftBC, rightBC;
  const ordinal_type interiorStencilLength, cornerStencilLength;
  ordinal_type numInterior;
  size_type numEntries;

  // Matrix views
  row_map_view_t rowmap;
  cols_view_t columns;
  scalar_view_t values;

  fill_1D_matrix_functor(const ordinal_type numNodes_, const int leftBC_, const int rightBC_,
                         const row_map_view_t rowmap_, const cols_view_t columns_, const scalar_view_t values_)
      : numNodes(numNodes_),
        leftBC(leftBC_),
        rightBC(rightBC_),
        interiorStencilLength(3),
        cornerStencilLength(2),
        rowmap(rowmap_),
        columns(columns_),
        values(values_) {
    if (numNodes == 1) {
      std::ostringstream os;
      os << "You need at least two points per direction to obtain a valid "
            "discretization !"
         << std::endl;
      throw std::runtime_error(os.str());
    }

    numInterior = numNodes - 2;
    numEntries  = numInterior * interiorStencilLength + 2 * cornerStencilLength;
  }

  void compute() {
    // Fill interior points
    if (0 < numInterior) {
      Kokkos::RangePolicy<execution_space, interiorTag> interiorPolicy(0, numInterior);
      Kokkos::parallel_for("Fill 1D matrix: interior points", interiorPolicy, *this);
    }

    // Fill exterior points a.k.a. boundary points
    Kokkos::RangePolicy<execution_space, exteriorTag> exteriorPolicy(0, 1);
    Kokkos::parallel_for("Fill 1D matrix: exterior points", exteriorPolicy, *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const interiorTag&, const ordinal_type idx) const {
    const ordinal_type rowIdx = idx + 1;  // Offset by one since first node has BC
    const size_type rowOffset = size_type(rowIdx) * interiorStencilLength + cornerStencilLength;

    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 3) = rowIdx - 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + 1;

    // Fill values
    values(rowOffset - 3) = -1.0;
    values(rowOffset - 2) = 2.0;
    values(rowOffset - 1) = -1.0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const exteriorTag&, const ordinal_type /*idx*/) const {
    // LeftBC
    rowmap(1) = 2;

    columns(0) = 0;
    columns(1) = 1;
    if (leftBC == 1) {
      values(0) = 1.0;
      values(1) = 0.0;
    } else {
      values(0) = 1.0;
      values(1) = -1.0;
    }

    // RightBC
    rowmap(numNodes) = numEntries;

    columns(numEntries - 2) = numNodes - 2;
    columns(numEntries - 1) = numNodes - 1;
    if (rightBC == 1) {
      values(numEntries - 2) = 0.0;
      values(numEntries - 1) = 1.0;
    } else {
      values(numEntries - 2) = -1.0;
      values(numEntries - 1) = 1.0;
    }
  }
};

template <typename CrsMatrix_t, typename mat_structure>
CrsMatrix_t generate_structured_matrix1D(const mat_structure& structure) {
  typedef typename CrsMatrix_t::StaticCrsGraphType graph_t;
  typedef typename CrsMatrix_t::row_map_type::non_const_type row_map_view_t;
  typedef typename CrsMatrix_t::index_type::non_const_type cols_view_t;
  typedef typename CrsMatrix_t::values_type::non_const_type scalar_view_t;
  typedef typename CrsMatrix_t::non_const_size_type size_type;
  typedef typename CrsMatrix_t::non_const_ordinal_type ordinal_type;

  // Extract geometric data
  const ordinal_type nx                    = structure(0, 0);
  const ordinal_type numNodes              = nx;
  const ordinal_type leftBC                = structure(0, 1);
  const ordinal_type rightBC               = structure(0, 2);
  const ordinal_type numInterior           = (nx - leftBC - rightBC);
  const ordinal_type numCorner             = leftBC + rightBC;
  const ordinal_type interiorStencilLength = 3, cornerStencilLength = 2;
  const size_type numEntries = numInterior * interiorStencilLength + numCorner * cornerStencilLength;

  // Create matrix data
  row_map_view_t rowmap_view("rowmap_view", numNodes + 1);
  cols_view_t columns_view("colsmap_view", numEntries);
  scalar_view_t values_view("values_view", numEntries);

  fill_1D_matrix_functor<CrsMatrix_t> fill_matrix(numNodes, leftBC, rightBC, rowmap_view, columns_view, values_view);
  fill_matrix.compute();

  graph_t static_graph(columns_view, rowmap_view);
  std::string name = "CrsMatrixFE";

  return CrsMatrix_t(name, numNodes, values_view, static_graph);

}  // generate_structured_matrix1D

template <class CrsMatrix_t>
struct fill_2D_matrix_functor {
  // Define types used by the CrsMatrix
  using execution_space = typename CrsMatrix_t::execution_space;
  using row_map_view_t  = typename CrsMatrix_t::row_map_type::non_const_type;
  using cols_view_t     = typename CrsMatrix_t::index_type::non_const_type;
  using scalar_view_t   = typename CrsMatrix_t::values_type::non_const_type;
  using size_type       = typename CrsMatrix_t::non_const_size_type;
  using ordinal_type    = typename CrsMatrix_t::non_const_ordinal_type;

  // Finite difference dispatch tags
  struct interiorFDTag {};

  struct xEdgeFDTag {};
  struct yEdgeFDTag {};

  struct cornerFDTag {};

  // Finite element dispatch tags
  struct interiorFETag {};

  struct xEdgeFETag {};
  struct yEdgeFETag {};

  struct cornerFETag {};

  // Internal variables and temporaries
  const int stencil_type;
  const ordinal_type nx, ny;
  const int leftBC, rightBC, bottomBC, topBC;

  // Matrix views
  row_map_view_t rowmap;
  cols_view_t columns;
  scalar_view_t values;

  ordinal_type interiorStencilLength, edgeStencilLength, cornerStencilLength;
  ordinal_type numInterior;
  ordinal_type numXEdge;
  ordinal_type numYEdge;
  ordinal_type numCorner;
  ordinal_type numEntriesPerGridRow;
  ordinal_type numEntriesBottomRow;
  size_type numEntries;

  fill_2D_matrix_functor(const int stencil_type_, const ordinal_type nx_, const ordinal_type ny_, const int leftBC_,
                         const int rightBC_, const int bottomBC_, const int topBC_, const row_map_view_t rowmap_,
                         const cols_view_t columns_, const scalar_view_t values_)
      : stencil_type(stencil_type_),
        nx(nx_),
        ny(ny_),
        leftBC(leftBC_),
        rightBC(rightBC_),
        bottomBC(bottomBC_),
        topBC(topBC_),
        rowmap(rowmap_),
        columns(columns_),
        values(values_) {
    if (nx == 1 || ny == 1) {
      std::ostringstream os;
      os << "You need at least two points per direction to obtain a valid "
            "discretization!"
         << std::endl;
      throw std::runtime_error(os.str());
    }

    if (stencil_type == FD) {
      interiorStencilLength = 5;
      edgeStencilLength     = 4;
      cornerStencilLength   = 3;
    } else if (stencil_type == FE) {
      interiorStencilLength = 9;
      edgeStencilLength     = 6;
      cornerStencilLength   = 4;
    }

    numInterior = (nx - 2) * (ny - 2);
    numXEdge    = nx - 2;
    numYEdge    = ny - 2;
    numCorner   = 4;

    numEntriesPerGridRow = (nx - 2) * interiorStencilLength + 2 * edgeStencilLength;

    numEntriesBottomRow = (nx - 2) * edgeStencilLength + 2 * cornerStencilLength;

    numEntries = numInterior * interiorStencilLength + (2 * numXEdge + 2 * numYEdge) * edgeStencilLength +
                 numCorner * cornerStencilLength;
  }

  void compute() {
    // Fill interior points
    if (0 < numInterior) {
      if (stencil_type == FD) {
        Kokkos::RangePolicy<execution_space, interiorFDTag> policy(0, numInterior);
        Kokkos::parallel_for("Fill 2D FD matrix: interior points", policy, *this);
      } else if (stencil_type == FE) {
        Kokkos::RangePolicy<execution_space, interiorFETag> policy(0, numInterior);
        Kokkos::parallel_for("Fill 2D FE matrix: interior points", policy, *this);
      }
    }

    // Fill x-edge points
    if (0 < numXEdge) {
      if (stencil_type == FD) {
        Kokkos::RangePolicy<execution_space, xEdgeFDTag> policy(0, numXEdge);
        Kokkos::parallel_for("Fill 2D FD matrix: x-edge points", policy, *this);
      } else if (stencil_type == FE) {
        Kokkos::RangePolicy<execution_space, xEdgeFETag> policy(0, numXEdge);
        Kokkos::parallel_for("Fill 2D FE matrix: x-edge points", policy, *this);
      }
    }

    // Fill y-edge points
    if (0 < numYEdge) {
      if (stencil_type == FD) {
        Kokkos::RangePolicy<execution_space, yEdgeFDTag> policy(0, numYEdge);
        Kokkos::parallel_for("Fill 2D FD matrix: y-edge points", policy, *this);
      } else if (stencil_type == FE) {
        Kokkos::RangePolicy<execution_space, yEdgeFETag> policy(0, numYEdge);
        Kokkos::parallel_for("Fill 2D FE matrix: y-edge points", policy, *this);
      }
    }

    // Fill corner points
    if (0 < numCorner) {
      if (stencil_type == FD) {
        Kokkos::RangePolicy<execution_space, cornerFDTag> policy(0, 1);
        Kokkos::parallel_for("Fill 2D FD matrix: corner points", policy, *this);
      } else if (stencil_type == FE) {
        Kokkos::RangePolicy<execution_space, cornerFETag> policy(0, 1);
        Kokkos::parallel_for("Fill 2D FE matrix: corner points", policy, *this);
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const interiorFDTag&, const ordinal_type idx) const {
    ordinal_type i, j;

    // Compute row index
    j                         = idx / (nx - 2);
    i                         = idx % (nx - 2);
    const ordinal_type rowIdx = (j + 1) * nx + i + 1;

    // Compute rowOffset
    const size_type rowOffset = size_type(j) * numEntriesPerGridRow + numEntriesBottomRow +
                                size_type(i + 1) * interiorStencilLength + edgeStencilLength;

    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - nx;
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + 1;
    columns(rowOffset - 1) = rowIdx + nx;

    // Fill values
    values(rowOffset - 5) = -1.0;
    values(rowOffset - 4) = -1.0;
    values(rowOffset - 3) = 4.0;
    values(rowOffset - 2) = -1.0;
    values(rowOffset - 1) = -1.0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const xEdgeFDTag&, const ordinal_type idx) const {
    /***************/
    /* Bottom edge */
    /***************/
    ordinal_type rowIdx = idx + 1;
    size_type rowOffset = size_type(idx + 1) * edgeStencilLength + cornerStencilLength;

    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + 1;
    columns(rowOffset - 1) = rowIdx + nx;
    if (bottomBC == 1) {
      // Fill values
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 3.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    /************/
    /* Top edge */
    /************/
    rowIdx    = (ny - 1) * nx + idx + 1;
    rowOffset = size_type(ny - 2) * numEntriesPerGridRow + numEntriesBottomRow +
                size_type(idx + 1) * edgeStencilLength + cornerStencilLength;

    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx - nx;
    columns(rowOffset - 3) = rowIdx - 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + 1;
    if (topBC == 1) {
      // Fill values
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 3.0;
      values(rowOffset - 1) = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const yEdgeFDTag&, const ordinal_type idx) const {
    /*************/
    /* Left edge */
    /*************/
    ordinal_type rowIdx = (idx + 1) * nx;
    size_type rowOffset = size_type(idx) * numEntriesPerGridRow + numEntriesBottomRow + edgeStencilLength;

    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx - nx;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + 1;
    columns(rowOffset - 1) = rowIdx + nx;
    if (leftBC == 1) {
      // Fill values
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 3.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    /**************/
    /* Right edge */
    /**************/
    rowIdx             = (idx + 2) * nx - 1;
    rowOffset          = size_type(idx + 1) * numEntriesPerGridRow + numEntriesBottomRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx - nx;
    columns(rowOffset - 3) = rowIdx - 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + nx;
    if (rightBC == 1) {
      // Fill values
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 3.0;
      values(rowOffset - 1) = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const cornerFDTag&, const ordinal_type /*idx*/) const {
    // Bottom-left corner
    ordinal_type rowIdx = 0;
    size_type rowOffset = cornerStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + 1;
    columns(rowOffset - 1) = rowIdx + nx;
    if (bottomBC == 1 || leftBC == 1) {
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      values(rowOffset - 3) = 2.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    // Bottom-right corner
    rowIdx             = nx - 1;
    rowOffset          = numEntriesBottomRow;
    rowmap(rowIdx + 1) = rowOffset;

    columns(rowOffset - 3) = rowIdx - 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + nx;
    if (bottomBC == 1 || rightBC == 1) {
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 2.0;
      values(rowOffset - 1) = -1.0;
    }

    // Top-left corner
    rowIdx             = (ny - 1) * nx;
    rowOffset          = size_type(ny - 2) * numEntriesPerGridRow + numEntriesBottomRow + cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    columns(rowOffset - 3) = rowIdx - nx;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + 1;
    if (topBC == 1 || leftBC == 1) {
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 2.0;
      values(rowOffset - 1) = -1.0;
    }

    // Top-right corner
    rowIdx             = ny * nx - 1;
    rowOffset          = numEntries;
    rowmap(rowIdx + 1) = rowOffset;

    columns(rowOffset - 3) = rowIdx - nx;
    columns(rowOffset - 2) = rowIdx - 1;
    columns(rowOffset - 1) = rowIdx;
    if (topBC == 1 || rightBC == 1) {
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 1.0;
    } else {
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = 2.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const interiorFETag&, const ordinal_type idx) const {
    ordinal_type i, j;

    // Compute row index
    j                         = idx / (nx - 2);
    i                         = idx % (nx - 2);
    const ordinal_type rowIdx = (j + 1) * nx + i + 1;

    // Compute rowOffset
    const size_type rowOffset = size_type(j) * numEntriesPerGridRow + numEntriesBottomRow +
                                size_type(i + 1) * interiorStencilLength + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 9) = rowIdx - nx - 1;
    columns(rowOffset - 8) = rowIdx - nx;
    columns(rowOffset - 7) = rowIdx - nx + 1;
    columns(rowOffset - 6) = rowIdx - 1;
    columns(rowOffset - 5) = rowIdx;
    columns(rowOffset - 4) = rowIdx + 1;
    columns(rowOffset - 3) = rowIdx + nx - 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + nx + 1;

    // Fill values
    values(rowOffset - 9) = -2.0;
    values(rowOffset - 8) = -2.0;
    values(rowOffset - 7) = -2.0;
    values(rowOffset - 6) = -2.0;
    values(rowOffset - 5) = 16.0;
    values(rowOffset - 4) = -2.0;
    values(rowOffset - 3) = -2.0;
    values(rowOffset - 2) = -2.0;
    values(rowOffset - 1) = -2.0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const xEdgeFETag&, const ordinal_type idx) const {
    /***************/
    /* Bottom edge */
    /***************/
    ordinal_type rowIdx = idx + 1;
    size_type rowOffset = size_type(idx + 1) * edgeStencilLength + cornerStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    // Fill column indices
    columns(rowOffset - 6) = rowIdx - 1;
    columns(rowOffset - 5) = rowIdx;
    columns(rowOffset - 4) = rowIdx + 1;
    columns(rowOffset - 3) = rowIdx + nx - 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + nx + 1;
    if (bottomBC == 1) {
      // Fill values
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 1.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = 8.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -2.0;
      values(rowOffset - 2) = -2.0;
      values(rowOffset - 1) = -2.0;
    }

    /************/
    /* Top edge */
    /************/
    rowIdx    = (ny - 1) * nx + idx + 1;
    rowOffset = size_type(ny - 2) * numEntriesPerGridRow + numEntriesBottomRow +
                size_type(idx + 1) * edgeStencilLength + cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 6) = rowIdx - nx - 1;
    columns(rowOffset - 5) = rowIdx - nx;
    columns(rowOffset - 4) = rowIdx - nx + 1;
    columns(rowOffset - 3) = rowIdx - 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + 1;
    if (topBC == 1) {
      // Fill values
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 6) = -2.0;
      values(rowOffset - 5) = -2.0;
      values(rowOffset - 4) = -2.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 8.0;
      values(rowOffset - 1) = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const yEdgeFETag&, const ordinal_type idx) const {
    /*************/
    /* Left edge */
    /*************/
    ordinal_type rowIdx = (idx + 1) * nx;
    size_type rowOffset = size_type(idx) * numEntriesPerGridRow + numEntriesBottomRow + edgeStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    // Fill column indices
    columns(rowOffset - 6) = rowIdx - nx;
    columns(rowOffset - 5) = rowIdx - nx + 1;
    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + nx + 1;
    if (leftBC == 1) {
      // Fill values
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = -2.0;
      values(rowOffset - 4) = 8.0;
      values(rowOffset - 3) = -2.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -2.0;
    }

    /**************/
    /* Right edge */
    /**************/
    rowIdx             = (idx + 2) * nx - 1;
    rowOffset          = size_type(idx + 1) * numEntriesPerGridRow + numEntriesBottomRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 6) = rowIdx - nx - 1;
    columns(rowOffset - 5) = rowIdx - nx;
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + nx - 1;
    columns(rowOffset - 1) = rowIdx + nx;
    if (rightBC == 1) {
      // Fill values
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 6) = -2.0;
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -2.0;
      values(rowOffset - 3) = 8.0;
      values(rowOffset - 2) = -2.0;
      values(rowOffset - 1) = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const cornerFETag&, const ordinal_type /*idx*/) const {
    // Bottom-left corner
    ordinal_type rowIdx = 0;
    size_type rowOffset = cornerStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + nx + 1;
    if (bottomBC == 1 || leftBC == 1) {
      values(rowOffset - 4) = 1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      values(rowOffset - 4) = 4.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -2.0;
    }

    // Bottom-right corner
    rowIdx             = nx - 1;
    rowOffset          = numEntriesBottomRow;
    rowmap(rowIdx + 1) = rowOffset;

    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + nx - 1;
    columns(rowOffset - 1) = rowIdx + nx;
    if (bottomBC == 1 || rightBC == 1) {
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 4.0;
      values(rowOffset - 2) = -2.0;
      values(rowOffset - 1) = -1.0;
    }

    // Top-left corner
    rowIdx             = (ny - 1) * nx;
    rowOffset          = size_type(ny - 2) * numEntriesPerGridRow + numEntriesBottomRow + cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    columns(rowOffset - 4) = rowIdx - nx;
    columns(rowOffset - 3) = rowIdx - nx + 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + 1;
    if (topBC == 1 || leftBC == 1) {
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -2.0;
      values(rowOffset - 2) = 4.0;
      values(rowOffset - 1) = -1.0;
    }

    // Top-right corner
    rowIdx             = ny * nx - 1;
    rowOffset          = numEntries;
    rowmap(rowIdx + 1) = rowOffset;

    columns(rowOffset - 4) = rowIdx - nx - 1;
    columns(rowOffset - 3) = rowIdx - nx;
    columns(rowOffset - 2) = rowIdx - 1;
    columns(rowOffset - 1) = rowIdx;
    if (topBC == 1 || rightBC == 1) {
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 1.0;
    } else {
      values(rowOffset - 4) = -2.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = 4.0;
    }
  }
};

template <typename CrsMatrix_t, typename mat_structure>
CrsMatrix_t generate_structured_matrix2D(const std::string stencil, const mat_structure& structure) {
  typedef typename CrsMatrix_t::StaticCrsGraphType graph_t;
  typedef typename CrsMatrix_t::row_map_type::non_const_type row_map_view_t;
  typedef typename CrsMatrix_t::index_type::non_const_type cols_view_t;
  typedef typename CrsMatrix_t::values_type::non_const_type scalar_view_t;
  typedef typename CrsMatrix_t::non_const_size_type size_type;
  typedef typename CrsMatrix_t::non_const_ordinal_type ordinal_type;

  int stencil_type = 0;
  if (stencil == "FD") {
    stencil_type = FD;
  } else if (stencil == "FE") {
    stencil_type = FE;
  } else {
    std::ostringstream os;
    os << "Test::generate_structured_matrix2D only accepts stencil: FD and "
          "FEM, you passed: "
       << stencil << " !" << std::endl;
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Extract geometric data
  const ordinal_type nx              = structure(0, 0);
  const ordinal_type ny              = structure(1, 0);
  const ordinal_type numNodes        = ny * nx;
  const ordinal_type leftBC          = structure(0, 1);
  const ordinal_type rightBC         = structure(0, 2);
  const ordinal_type bottomBC        = structure(1, 1);
  const ordinal_type topBC           = structure(1, 2);
  const ordinal_type numInterior     = (nx - 2) * (ny - 2);
  const ordinal_type numEdge         = 2 * (nx - 2) + 2 * (ny - 2);
  const ordinal_type numCorner       = 4;
  ordinal_type interiorStencilLength = 0, edgeStencilLength = 0, cornerStencilLength = 0;

  if (stencil_type == FD) {
    interiorStencilLength = 5;
    edgeStencilLength     = 4;
    cornerStencilLength   = 3;
  } else if (stencil_type == FE) {
    interiorStencilLength = 9;
    edgeStencilLength     = 6;
    cornerStencilLength   = 4;
  }

  const size_type numEntries =
      numInterior * interiorStencilLength + numEdge * edgeStencilLength + numCorner * cornerStencilLength;

  // Create matrix data
  row_map_view_t rowmap_view("rowmap_view", numNodes + 1);
  cols_view_t columns_view("colsmap_view", numEntries);
  scalar_view_t values_view("values_view", numEntries);

  fill_2D_matrix_functor<CrsMatrix_t> fill_2D_matrix(stencil_type, nx, ny, leftBC, rightBC, bottomBC, topBC,
                                                     rowmap_view, columns_view, values_view);

  fill_2D_matrix.compute();

  graph_t static_graph(columns_view, rowmap_view);
  std::string name;
  if (stencil_type == FD) {
    name = "CrsMatrixFD";
  } else if (stencil_type == FE) {
    name = "CrsMatrixFE";
  }

  return CrsMatrix_t(name, numNodes, values_view, static_graph);

}  // generate_structured_matrix2D

template <class CrsMatrix_t>
struct fill_3D_matrix_functor {
  // Define types used by the CrsMatrix
  typedef typename CrsMatrix_t::execution_space execution_space;
  typedef typename CrsMatrix_t::row_map_type::non_const_type row_map_view_t;
  typedef typename CrsMatrix_t::index_type::non_const_type cols_view_t;
  typedef typename CrsMatrix_t::values_type::non_const_type scalar_view_t;
  typedef typename CrsMatrix_t::non_const_size_type size_type;
  typedef typename CrsMatrix_t::non_const_ordinal_type ordinal_type;

  // Finite Difference dispatch tags
  struct interiorFDTag {};

  struct xFaceFDTag {};
  struct yFaceFDTag {};
  struct zFaceFDTag {};

  struct xEdgeFDTag {};
  struct yEdgeFDTag {};
  struct zEdgeFDTag {};

  struct cornerFDTag {};

  // Finite Element dispatch tags
  struct interiorFETag {};

  struct xFaceFETag {};
  struct yFaceFETag {};
  struct zFaceFETag {};

  struct xEdgeFETag {};
  struct yEdgeFETag {};
  struct zEdgeFETag {};

  struct cornerFETag {};

  // Internal variables and temporaries
  const int stencil_type;
  const ordinal_type nx, ny, nz;
  const int leftBC, rightBC, frontBC, backBC, bottomBC, topBC;

  // Matrix views
  row_map_view_t rowmap;
  cols_view_t columns;
  scalar_view_t values;

  size_type numEntries;
  ordinal_type numInterior;
  ordinal_type numXFace;
  ordinal_type numYFace;
  ordinal_type numZFace;
  ordinal_type numXEdge;
  ordinal_type numYEdge;
  ordinal_type numZEdge;

  ordinal_type interiorStencilLength;
  ordinal_type faceStencilLength;
  ordinal_type edgeStencilLength;
  ordinal_type cornerStencilLength;
  ordinal_type numEntriesPerGridPlane;
  ordinal_type numEntriesBottomPlane;
  ordinal_type numEntriesPerGridRow;
  ordinal_type numEntriesFrontRow;
  ordinal_type numEntriesBottomFrontRow;

  fill_3D_matrix_functor(const int stencil_type_, const ordinal_type nx_, const ordinal_type ny_,
                         const ordinal_type nz_, const int leftBC_, const int rightBC_, const int frontBC_,
                         const int backBC_, const int bottomBC_, const int topBC_, const row_map_view_t rowmap_,
                         const cols_view_t columns_, const scalar_view_t values_)
      : stencil_type(stencil_type_),
        nx(nx_),
        ny(ny_),
        nz(nz_),
        leftBC(leftBC_),
        rightBC(rightBC_),
        frontBC(frontBC_),
        backBC(backBC_),
        bottomBC(bottomBC_),
        topBC(topBC_),
        rowmap(rowmap_),
        columns(columns_),
        values(values_) {
    if (stencil_type == FD) {
      interiorStencilLength = 7;
      faceStencilLength     = 6;
      edgeStencilLength     = 5;
      cornerStencilLength   = 4;
    } else if (stencil_type == FE) {
      interiorStencilLength = 27;
      faceStencilLength     = 18;
      edgeStencilLength     = 12;
      cornerStencilLength   = 8;
    }

    numInterior = (nx - 2) * (ny - 2) * (nz - 2);
    numXFace    = (ny - 2) * (nz - 2);
    numYFace    = (nx - 2) * (nz - 2);
    numZFace    = (nx - 2) * (ny - 2);
    numXEdge    = nx - 2;
    numYEdge    = ny - 2;
    numZEdge    = nz - 2;

    numEntries = numInterior * interiorStencilLength + 2 * (numXFace + numYFace + numZFace) * faceStencilLength +
                 4 * (numXEdge + numYEdge + numZEdge) * edgeStencilLength + 8 * cornerStencilLength;
    numEntriesPerGridPlane = numZFace * interiorStencilLength + 2 * numXEdge * faceStencilLength +
                             2 * numYEdge * faceStencilLength + 4 * edgeStencilLength;
    ;
    numEntriesBottomPlane = numZFace * faceStencilLength + 2 * numXEdge * edgeStencilLength +
                            2 * numYEdge * edgeStencilLength + 4 * cornerStencilLength;
    ;
    numEntriesPerGridRow     = numXEdge * interiorStencilLength + 2 * faceStencilLength;
    numEntriesFrontRow       = numXEdge * faceStencilLength + 2 * edgeStencilLength;
    numEntriesBottomFrontRow = numXEdge * edgeStencilLength + 2 * cornerStencilLength;
  }

  void compute() {
    // Fill interior points
    if (0 < numInterior) {
      if (stencil_type == FD) {
        Kokkos::RangePolicy<execution_space, interiorFDTag> policy(0, numInterior);
        Kokkos::parallel_for("Fill 3D FD matrix: interior points", policy, *this);
      } else if (stencil_type == FE) {
        Kokkos::RangePolicy<execution_space, interiorFETag> policy(0, numInterior);
        Kokkos::parallel_for("Fill 3D FE matrix: interior points", policy, *this);
      }
    }

    // Fill x-faces
    if (0 < numXFace) {
      if (stencil_type == FD) {
        Kokkos::RangePolicy<execution_space, xFaceFDTag> policy(0, numXFace);
        Kokkos::parallel_for("Fill 3D FD matrix: x-face points", policy, *this);
      } else if (stencil_type == FE) {
        Kokkos::RangePolicy<execution_space, xFaceFETag> policy(0, numXFace);
        Kokkos::parallel_for("Fill 3D FE matrix: x-face points", policy, *this);
      }
    }

    // Fill y-faces
    if (0 < numYFace) {
      if (stencil_type == FD) {
        Kokkos::RangePolicy<execution_space, yFaceFDTag> policy(0, numYFace);
        Kokkos::parallel_for("Fill 3D FD matrix: y-face points", policy, *this);
      } else if (stencil_type == FE) {
        Kokkos::RangePolicy<execution_space, yFaceFETag> policy(0, numYFace);
        Kokkos::parallel_for("Fill 3D FE matrix: y-face points", policy, *this);
      }
    }

    // Fill z-faces
    if (0 < numZFace) {
      if (stencil_type == FD) {
        Kokkos::RangePolicy<execution_space, zFaceFDTag> policy(0, numZFace);
        Kokkos::parallel_for("Fill 3D FD matrix: z-face points", policy, *this);
      } else if (stencil_type == FE) {
        Kokkos::RangePolicy<execution_space, zFaceFETag> policy(0, numZFace);
        Kokkos::parallel_for("Fill 3D FE matrix: z-face points", policy, *this);
      }
    }

    // Fill x-edges
    if (0 < numXEdge) {
      if (stencil_type == FD) {
        Kokkos::RangePolicy<execution_space, xEdgeFDTag> policy(0, numXEdge);
        Kokkos::parallel_for("Fill 3D FD matrix: x-edge points", policy, *this);
      } else if (stencil_type == FE) {
        Kokkos::RangePolicy<execution_space, xEdgeFETag> policy(0, numXEdge);
        Kokkos::parallel_for("Fill 3D FE matrix: x-edge points", policy, *this);
      }
    }

    // Fill y-edges
    if (0 < numYEdge) {
      if (stencil_type == FD) {
        Kokkos::RangePolicy<execution_space, yEdgeFDTag> policy(0, numYEdge);
        Kokkos::parallel_for("Fill 3D FD matrix: y-edge points", policy, *this);
      } else if (stencil_type == FE) {
        Kokkos::RangePolicy<execution_space, yEdgeFETag> policy(0, numYEdge);
        Kokkos::parallel_for("Fill 3D FE matrix: y-edge points", policy, *this);
      }
    }

    // Fill z-edges
    if (0 < numZEdge) {
      if (stencil_type == FD) {
        Kokkos::RangePolicy<execution_space, zEdgeFDTag> policy(0, numZEdge);
        Kokkos::parallel_for("Fill 3D FD matrix: z-edge points", policy, *this);
      } else if (stencil_type == FE) {
        Kokkos::RangePolicy<execution_space, zEdgeFETag> policy(0, numZEdge);
        Kokkos::parallel_for("Fill 3D FE matrix: z-edge points", policy, *this);
      }
    }

    if (stencil_type == FD) {
      Kokkos::RangePolicy<execution_space, cornerFDTag> policy(0, 1);
      Kokkos::parallel_for("Fill 3D FD matrix: corner points", policy, *this);
    } else if (stencil_type == FE) {
      Kokkos::RangePolicy<execution_space, cornerFETag> policy(0, 1);
      Kokkos::parallel_for("Fill 3D FE matrix: corner points", policy, *this);
    }
  }  // compute()

  KOKKOS_INLINE_FUNCTION
  void operator()(const interiorFDTag&, const ordinal_type idx) const {
    // Compute row index
    const ordinal_type k      = idx / ((ny - 2) * (nx - 2));
    const ordinal_type rem    = idx % ((ny - 2) * (nx - 2));
    const ordinal_type j      = rem / (nx - 2);
    const ordinal_type i      = rem % (nx - 2);
    const ordinal_type rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i + 1;

    // Compute rowOffset
    const size_type rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane +
                                size_type(j) * numEntriesPerGridRow + numEntriesFrontRow +
                                size_type(i + 1) * interiorStencilLength + faceStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 7) = rowIdx - ny * nx;
    columns(rowOffset - 6) = rowIdx - nx;
    columns(rowOffset - 5) = rowIdx - 1;
    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;

    // Fill values
    values(rowOffset - 7) = -1.0;
    values(rowOffset - 6) = -1.0;
    values(rowOffset - 5) = -1.0;
    values(rowOffset - 4) = 6.0;
    values(rowOffset - 3) = -1.0;
    values(rowOffset - 2) = -1.0;
    values(rowOffset - 1) = -1.0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const xFaceFDTag&, const ordinal_type idx) const {
    /*******************/
    /*   x == 0 face   */
    /*******************/
    // Compute row index
    ordinal_type k      = idx / (ny - 2);
    ordinal_type j      = idx % (ny - 2);
    ordinal_type i      = 0;
    ordinal_type rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i;

    // Compute rowOffset
    size_type rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane +
                          size_type(j) * numEntriesPerGridRow + numEntriesFrontRow + faceStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 6) = rowIdx - ny * nx;
    columns(rowOffset - 5) = rowIdx - nx;
    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (leftBC == 1) {
      // Fill values
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = 5.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    /********************/
    /*   x == nx face   */
    /********************/
    // Compute row index
    k      = idx / (ny - 2);
    j      = idx % (ny - 2);
    i      = nx - 1;
    rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane +
                size_type(j + 1) * numEntriesPerGridRow + numEntriesFrontRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 6) = rowIdx - ny * nx;
    columns(rowOffset - 5) = rowIdx - nx;
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (rightBC == 1) {
      // Fill values
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 5.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const yFaceFDTag&, const ordinal_type idx) const {
    /*******************/
    /*   y == 0 face   */
    /*******************/
    // Compute row index
    ordinal_type k      = idx / (nx - 2);
    ordinal_type i      = idx % (nx - 2);
    ordinal_type rowIdx = (k + 1) * ny * nx + i + 1;

    // Compute rowOffset
    size_type rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane +
                          size_type(i + 1) * faceStencilLength + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 6) = rowIdx - ny * nx;
    columns(rowOffset - 5) = rowIdx - 1;
    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (frontBC == 1) {
      // Fill values
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = 5.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    /********************/
    /*   y == ny face   */
    /********************/
    // Compute row index
    ordinal_type j = ny - 2;
    k              = idx / (nx - 2);
    i              = idx % (nx - 2);
    rowIdx         = (k + 1) * ny * nx + (j + 1) * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesPerGridRow +
                numEntriesFrontRow + size_type(i + 1) * faceStencilLength + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 6) = rowIdx - ny * nx;
    columns(rowOffset - 5) = rowIdx - nx;
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + 1;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (backBC == 1) {
      // Fill values
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 5.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const zFaceFDTag&, const ordinal_type idx) const {
    /*******************/
    /*   z == 0 face   */
    /*******************/
    // Compute row index
    ordinal_type j      = idx / (nx - 2);
    ordinal_type i      = idx % (nx - 2);
    ordinal_type rowIdx = (j + 1) * nx + i + 1;

    // Compute rowOffset
    size_type rowOffset = size_type(j) * numEntriesFrontRow + numEntriesBottomFrontRow +
                          size_type(i + 1) * faceStencilLength + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 6) = rowIdx - nx;
    columns(rowOffset - 5) = rowIdx - 1;
    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (bottomBC == 1) {
      // Fill values
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = 5.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    /********************/
    /*   z == nz face   */
    /********************/
    // Compute row index
    ordinal_type k = nz - 2;
    j              = idx / (nx - 2);
    i              = idx % (nx - 2);
    rowIdx         = (k + 1) * ny * nx + (j + 1) * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesFrontRow +
                numEntriesBottomFrontRow + size_type(i + 1) * faceStencilLength + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 6) = rowIdx - ny * nx;
    columns(rowOffset - 5) = rowIdx - nx;
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + 1;
    columns(rowOffset - 1) = rowIdx + nx;
    if (topBC == 1) {
      // Fill values
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 5.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const xEdgeFDTag&, const ordinal_type idx) const {
    // Compute row index
    ordinal_type i      = idx;
    ordinal_type rowIdx = i + 1;

    // Compute rowOffset
    size_type rowOffset = size_type(i + 1) * edgeStencilLength + cornerStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - 1;
    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (bottomBC == 1 || frontBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = 4.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    // Compute row index
    ordinal_type j = ny - 2;
    i              = idx;
    rowIdx         = (j + 1) * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(j) * numEntriesFrontRow + numEntriesBottomFrontRow + size_type(i + 1) * edgeStencilLength +
                cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - 1;
    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (bottomBC == 1 || backBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = 4.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    // Compute row index
    ordinal_type k = nz - 2;
    i              = idx;
    rowIdx         = (k + 1) * ny * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(i + 1) * edgeStencilLength +
                cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - ny * nx;
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + 1;
    columns(rowOffset - 1) = rowIdx + nx;
    if (topBC == 1 || frontBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 4.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    // Compute row index
    k      = nz - 2;
    j      = ny - 2;
    i      = idx;
    rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesFrontRow +
                numEntriesBottomFrontRow + size_type(i + 1) * edgeStencilLength + cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - ny * nx;
    columns(rowOffset - 4) = rowIdx - nx;
    columns(rowOffset - 3) = rowIdx - 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + 1;
    if (topBC == 1 || backBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 4.0;
      values(rowOffset - 1) = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const yEdgeFDTag&, const ordinal_type idx) const {
    // Compute row index
    ordinal_type j      = idx;
    ordinal_type rowIdx = (j + 1) * nx;

    // Compute rowOffset
    size_type rowOffset = size_type(j) * numEntriesFrontRow + numEntriesBottomFrontRow + edgeStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - nx;
    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (bottomBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = 4.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    // Compute row index
    j              = idx;
    ordinal_type i = nx - 1;
    rowIdx         = (j + 1) * nx + i;

    // Compute rowOffset
    rowOffset          = size_type(j + 1) * numEntriesFrontRow + numEntriesBottomFrontRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - nx;
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (bottomBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 4.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    // Compute row index
    ordinal_type k = nz - 2;
    j              = idx;
    rowIdx         = (k + 1) * ny * nx + (j + 1) * nx;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesFrontRow +
                numEntriesBottomFrontRow + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - ny * nx;
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + 1;
    columns(rowOffset - 1) = rowIdx + nx;
    if (topBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 4.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    // Compute row index
    k      = nz - 2;
    j      = idx;
    i      = nx - 1;
    rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j + 1) * numEntriesFrontRow +
                numEntriesBottomFrontRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - ny * nx;
    columns(rowOffset - 4) = rowIdx - nx;
    columns(rowOffset - 3) = rowIdx - 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + nx;
    if (topBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 4.0;
      values(rowOffset - 1) = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const zEdgeFDTag&, const ordinal_type idx) const {
    // Compute row index
    ordinal_type k      = idx;
    ordinal_type rowIdx = (k + 1) * ny * nx;

    // Compute rowOffset
    size_type rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + edgeStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - ny * nx;
    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (frontBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = 4.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    // Compute row index
    k              = idx;
    ordinal_type i = nx - 2;
    rowIdx         = (k + 1) * ny * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(i) * faceStencilLength +
                2 * edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - ny * nx;
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (frontBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 4.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    // Compute row index
    k              = idx;
    ordinal_type j = ny - 2;
    rowIdx         = (k + 1) * ny * nx + (j + 1) * nx;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesPerGridRow +
                numEntriesFrontRow + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - ny * nx;
    columns(rowOffset - 4) = rowIdx - nx;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + 1;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (backBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 4.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    // Compute row index
    k      = idx;
    j      = ny - 2;
    i      = nx - 2;
    rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesPerGridRow +
                numEntriesFrontRow + size_type(i) * faceStencilLength + 2 * edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 5) = rowIdx - ny * nx;
    columns(rowOffset - 4) = rowIdx - nx;
    columns(rowOffset - 3) = rowIdx - 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (backBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 4.0;
      values(rowOffset - 1) = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const cornerFDTag&, const ordinal_type /*idx*/) const {
    // Bottom corners
    ordinal_type rowIdx = 0;
    size_type rowOffset = cornerStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (bottomBC == 1 || frontBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 4) = 1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 4) = 3.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    rowIdx             = nx - 1;
    rowOffset          = numEntriesBottomFrontRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (bottomBC == 1 || frontBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 3.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    rowIdx             = (ny - 1) * nx;
    rowOffset          = size_type(ny - 2) * numEntriesFrontRow + numEntriesBottomFrontRow + cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx - nx;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + 1;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (bottomBC == 1 || backBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 3.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    rowIdx             = ny * nx - 1;
    rowOffset          = numEntriesBottomPlane;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx - nx;
    columns(rowOffset - 3) = rowIdx - 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (bottomBC == 1 || backBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 3.0;
      values(rowOffset - 1) = -1.0;
    }

    rowIdx             = (nz - 1) * ny * nx;
    rowOffset          = size_type(nz - 2) * numEntriesPerGridPlane + numEntriesBottomPlane + cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx - ny * nx;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + 1;
    columns(rowOffset - 1) = rowIdx + nx;
    if (topBC == 1 || frontBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 3.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    rowIdx             = (nz - 1) * ny * nx + nx - 1;
    rowOffset          = size_type(nz - 2) * numEntriesPerGridPlane + numEntriesBottomPlane + numEntriesBottomFrontRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx - ny * nx;
    columns(rowOffset - 3) = rowIdx - 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + nx;
    if (topBC == 1 || frontBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 3.0;
      values(rowOffset - 1) = -1.0;
    }

    rowIdx    = nz * ny * nx - nx;
    rowOffset = size_type(nz - 2) * numEntriesPerGridPlane + numEntriesBottomPlane + (ny - 2) * numEntriesFrontRow +
                numEntriesBottomFrontRow + cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx - ny * nx;
    columns(rowOffset - 3) = rowIdx - nx;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + 1;
    if (topBC == 1 || backBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 3.0;
      values(rowOffset - 1) = -1.0;
    }

    rowIdx             = nz * ny * nx - 1;
    rowOffset          = numEntries;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 4) = rowIdx - ny * nx;
    columns(rowOffset - 3) = rowIdx - nx;
    columns(rowOffset - 2) = rowIdx - 1;
    columns(rowOffset - 1) = rowIdx;
    if (topBC == 1 || backBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 1.0;
    } else {
      // Fill values
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = 3.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const interiorFETag&, const ordinal_type idx) const {
    // Compute row index
    const ordinal_type k      = idx / ((ny - 2) * (nx - 2));
    const ordinal_type rem    = idx % ((ny - 2) * (nx - 2));
    const ordinal_type j      = rem / (nx - 2);
    const ordinal_type i      = rem % (nx - 2);
    const ordinal_type rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i + 1;

    // Compute rowOffset
    const size_type rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane +
                                size_type(j) * numEntriesPerGridRow + numEntriesFrontRow +
                                size_type(i + 1) * interiorStencilLength + faceStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 27) = rowIdx - ny * nx - nx - 1;
    columns(rowOffset - 26) = rowIdx - ny * nx - nx;
    columns(rowOffset - 25) = rowIdx - ny * nx - nx + 1;
    columns(rowOffset - 24) = rowIdx - ny * nx - 1;
    columns(rowOffset - 23) = rowIdx - ny * nx;
    columns(rowOffset - 22) = rowIdx - ny * nx + 1;
    columns(rowOffset - 21) = rowIdx - ny * nx + nx - 1;
    columns(rowOffset - 20) = rowIdx - ny * nx + nx;
    columns(rowOffset - 19) = rowIdx - ny * nx + nx + 1;
    columns(rowOffset - 18) = rowIdx - nx - 1;
    columns(rowOffset - 17) = rowIdx - nx;
    columns(rowOffset - 16) = rowIdx - nx + 1;
    columns(rowOffset - 15) = rowIdx - 1;
    columns(rowOffset - 14) = rowIdx;
    columns(rowOffset - 13) = rowIdx + 1;
    columns(rowOffset - 12) = rowIdx + nx - 1;
    columns(rowOffset - 11) = rowIdx + nx;
    columns(rowOffset - 10) = rowIdx + nx + 1;
    columns(rowOffset - 9)  = rowIdx + nx * ny - nx - 1;
    columns(rowOffset - 8)  = rowIdx + nx * ny - nx;
    columns(rowOffset - 7)  = rowIdx + nx * ny - nx + 1;
    columns(rowOffset - 6)  = rowIdx + nx * ny - 1;
    columns(rowOffset - 5)  = rowIdx + nx * ny;
    columns(rowOffset - 4)  = rowIdx + nx * ny + 1;
    columns(rowOffset - 3)  = rowIdx + nx * ny + nx - 1;
    columns(rowOffset - 2)  = rowIdx + nx * ny + nx;
    columns(rowOffset - 1)  = rowIdx + nx * ny + nx + 1;

    // Fill values
    values(rowOffset - 27) = -1.0;
    values(rowOffset - 26) = -2.0;
    values(rowOffset - 25) = -1.0;
    values(rowOffset - 24) = -2.0;
    values(rowOffset - 23) = 0.0;
    values(rowOffset - 22) = -2.0;
    values(rowOffset - 21) = -1.0;
    values(rowOffset - 20) = -2.0;
    values(rowOffset - 19) = -1.0;
    values(rowOffset - 18) = -2.0;
    values(rowOffset - 17) = 0.0;
    values(rowOffset - 16) = -2.0;
    values(rowOffset - 15) = 0.0;
    values(rowOffset - 14) = 32.0;
    values(rowOffset - 13) = 0.0;
    values(rowOffset - 12) = -2.0;
    values(rowOffset - 11) = 0.0;
    values(rowOffset - 10) = -2.0;
    values(rowOffset - 9)  = -1.0;
    values(rowOffset - 8)  = -2.0;
    values(rowOffset - 7)  = -1.0;
    values(rowOffset - 6)  = -2.0;
    values(rowOffset - 5)  = 0.0;
    values(rowOffset - 4)  = -2.0;
    values(rowOffset - 3)  = -1.0;
    values(rowOffset - 2)  = -2.0;
    values(rowOffset - 1)  = -1.0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const xFaceFETag&, const ordinal_type idx) const {
    /*******************/
    /*   x == 0 face   */
    /*******************/
    // Compute row index
    ordinal_type k      = idx / (ny - 2);
    ordinal_type j      = idx % (ny - 2);
    ordinal_type i      = 0;
    ordinal_type rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i;

    // Compute rowOffset
    size_type rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane +
                          size_type(j) * numEntriesPerGridRow + numEntriesFrontRow + faceStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 18) = rowIdx - ny * nx - nx;
    columns(rowOffset - 17) = rowIdx - ny * nx - nx + 1;
    columns(rowOffset - 16) = rowIdx - ny * nx;
    columns(rowOffset - 15) = rowIdx - ny * nx + 1;
    columns(rowOffset - 14) = rowIdx - ny * nx + nx;
    columns(rowOffset - 13) = rowIdx - ny * nx + nx + 1;
    columns(rowOffset - 12) = rowIdx - nx;
    columns(rowOffset - 11) = rowIdx - nx + 1;
    columns(rowOffset - 10) = rowIdx;
    columns(rowOffset - 9)  = rowIdx + 1;
    columns(rowOffset - 8)  = rowIdx + nx;
    columns(rowOffset - 7)  = rowIdx + nx + 1;
    columns(rowOffset - 6)  = rowIdx + nx * ny - nx;
    columns(rowOffset - 5)  = rowIdx + nx * ny - nx + 1;
    columns(rowOffset - 4)  = rowIdx + nx * ny;
    columns(rowOffset - 3)  = rowIdx + nx * ny + 1;
    columns(rowOffset - 2)  = rowIdx + nx * ny + nx;
    columns(rowOffset - 1)  = rowIdx + nx * ny + nx + 1;
    if (leftBC == 1) {
      // Fill values
      values(rowOffset - 18) = 0.0;
      values(rowOffset - 17) = 0.0;
      values(rowOffset - 16) = 0.0;
      values(rowOffset - 15) = 0.0;
      values(rowOffset - 14) = 0.0;
      values(rowOffset - 13) = 0.0;
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 1.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 18) = -1.0;
      values(rowOffset - 17) = -1.0;
      values(rowOffset - 16) = 0.0;
      values(rowOffset - 15) = -2.0;
      values(rowOffset - 14) = -1.0;
      values(rowOffset - 13) = -1.0;
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = -2.0;
      values(rowOffset - 10) = 16.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = -2.0;
      values(rowOffset - 6)  = -1.0;
      values(rowOffset - 5)  = -1.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = -2.0;
      values(rowOffset - 2)  = -1.0;
      values(rowOffset - 1)  = -1.0;
    }

    /********************/
    /*   x == nx face   */
    /********************/
    // Compute row index
    k      = idx / (ny - 2);
    j      = idx % (ny - 2);
    i      = nx - 1;
    rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane +
                size_type(j + 1) * numEntriesPerGridRow + numEntriesFrontRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 18) = rowIdx - ny * nx - nx - 1;
    columns(rowOffset - 17) = rowIdx - ny * nx - nx;
    columns(rowOffset - 16) = rowIdx - ny * nx - 1;
    columns(rowOffset - 15) = rowIdx - ny * nx;
    columns(rowOffset - 14) = rowIdx - ny * nx + nx - 1;
    columns(rowOffset - 13) = rowIdx - ny * nx + nx;
    columns(rowOffset - 12) = rowIdx - nx - 1;
    columns(rowOffset - 11) = rowIdx - nx;
    columns(rowOffset - 10) = rowIdx - 1;
    columns(rowOffset - 9)  = rowIdx;
    columns(rowOffset - 8)  = rowIdx + nx - 1;
    columns(rowOffset - 7)  = rowIdx + nx;
    columns(rowOffset - 6)  = rowIdx + nx * ny - nx - 1;
    columns(rowOffset - 5)  = rowIdx + nx * ny - nx;
    columns(rowOffset - 4)  = rowIdx + nx * ny - 1;
    columns(rowOffset - 3)  = rowIdx + nx * ny;
    columns(rowOffset - 2)  = rowIdx + nx * ny + nx - 1;
    columns(rowOffset - 1)  = rowIdx + nx * ny + nx;
    if (rightBC == 1) {
      // Fill values
      values(rowOffset - 18) = 0.0;
      values(rowOffset - 17) = 0.0;
      values(rowOffset - 16) = 0.0;
      values(rowOffset - 15) = 0.0;
      values(rowOffset - 14) = 0.0;
      values(rowOffset - 13) = 0.0;
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 1.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 18) = -1.0;
      values(rowOffset - 17) = -1.0;
      values(rowOffset - 16) = -2.0;
      values(rowOffset - 15) = 0.0;
      values(rowOffset - 14) = -1.0;
      values(rowOffset - 13) = -1.0;
      values(rowOffset - 12) = -2.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 16.0;
      values(rowOffset - 8)  = -2.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = -1.0;
      values(rowOffset - 5)  = -1.0;
      values(rowOffset - 4)  = -2.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = -1.0;
      values(rowOffset - 1)  = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const yFaceFETag&, const ordinal_type idx) const {
    /*******************/
    /*   y == 0 face   */
    /*******************/
    // Compute row index
    ordinal_type k      = idx / (nx - 2);
    ordinal_type i      = idx % (nx - 2);
    ordinal_type rowIdx = (k + 1) * ny * nx + i + 1;

    // Compute rowOffset
    size_type rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane +
                          size_type(i + 1) * faceStencilLength + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 18) = rowIdx - ny * nx - 1;
    columns(rowOffset - 17) = rowIdx - ny * nx;
    columns(rowOffset - 16) = rowIdx - ny * nx + 1;
    columns(rowOffset - 15) = rowIdx - ny * nx + nx - 1;
    columns(rowOffset - 14) = rowIdx - ny * nx + nx;
    columns(rowOffset - 13) = rowIdx - ny * nx + nx + 1;
    columns(rowOffset - 12) = rowIdx - 1;
    columns(rowOffset - 11) = rowIdx;
    columns(rowOffset - 10) = rowIdx + 1;
    columns(rowOffset - 9)  = rowIdx + nx - 1;
    columns(rowOffset - 8)  = rowIdx + nx;
    columns(rowOffset - 7)  = rowIdx + nx + 1;
    columns(rowOffset - 6)  = rowIdx + nx * ny - 1;
    columns(rowOffset - 5)  = rowIdx + nx * ny;
    columns(rowOffset - 4)  = rowIdx + nx * ny + 1;
    columns(rowOffset - 3)  = rowIdx + nx * ny + nx - 1;
    columns(rowOffset - 2)  = rowIdx + nx * ny + nx;
    columns(rowOffset - 1)  = rowIdx + nx * ny + nx + 1;
    if (frontBC == 1) {
      // Fill values
      values(rowOffset - 18) = 0.0;
      values(rowOffset - 17) = 0.0;
      values(rowOffset - 16) = 0.0;
      values(rowOffset - 15) = 0.0;
      values(rowOffset - 14) = 0.0;
      values(rowOffset - 13) = 0.0;
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 1.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 18) = -1.0;
      values(rowOffset - 17) = 0.0;
      values(rowOffset - 16) = -1.0;
      values(rowOffset - 15) = -1.0;
      values(rowOffset - 14) = -2.0;
      values(rowOffset - 13) = -1.0;
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 16.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = -2.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = -2.0;
      values(rowOffset - 6)  = -1.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = -1.0;
      values(rowOffset - 3)  = -1.0;
      values(rowOffset - 2)  = -2.0;
      values(rowOffset - 1)  = -1.0;
    }

    /********************/
    /*   y == ny face   */
    /********************/
    // Compute row index
    ordinal_type j = ny - 2;
    k              = idx / (nx - 2);
    i              = idx % (nx - 2);
    rowIdx         = (k + 1) * ny * nx + (j + 1) * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesPerGridRow +
                numEntriesFrontRow + size_type(i + 1) * faceStencilLength + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 18) = rowIdx - ny * nx - nx - 1;
    columns(rowOffset - 17) = rowIdx - ny * nx - nx;
    columns(rowOffset - 16) = rowIdx - ny * nx - 1;
    columns(rowOffset - 15) = rowIdx - ny * nx - 1;
    columns(rowOffset - 14) = rowIdx - ny * nx;
    columns(rowOffset - 13) = rowIdx - ny * nx + 1;
    columns(rowOffset - 12) = rowIdx - nx - 1;
    columns(rowOffset - 11) = rowIdx - nx;
    columns(rowOffset - 10) = rowIdx - nx + 1;
    columns(rowOffset - 9)  = rowIdx - 1;
    columns(rowOffset - 8)  = rowIdx;
    columns(rowOffset - 7)  = rowIdx + 1;
    columns(rowOffset - 6)  = rowIdx + nx * ny - nx - 1;
    columns(rowOffset - 5)  = rowIdx + nx * ny - nx;
    columns(rowOffset - 4)  = rowIdx + nx * ny - nx + 1;
    columns(rowOffset - 3)  = rowIdx + nx * ny - 1;
    columns(rowOffset - 2)  = rowIdx + nx * ny;
    columns(rowOffset - 1)  = rowIdx + nx * ny + 1;
    if (backBC == 1) {
      // Fill values
      values(rowOffset - 18) = 0.0;
      values(rowOffset - 17) = 0.0;
      values(rowOffset - 16) = 0.0;
      values(rowOffset - 15) = 0.0;
      values(rowOffset - 14) = 0.0;
      values(rowOffset - 13) = 0.0;
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 1.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 18) = -1.0;
      values(rowOffset - 17) = -2.0;
      values(rowOffset - 16) = -1.0;
      values(rowOffset - 15) = -1.0;
      values(rowOffset - 14) = 0.0;
      values(rowOffset - 13) = -1.0;
      values(rowOffset - 12) = -2.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = -2.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 16.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = -1.0;
      values(rowOffset - 5)  = -2.0;
      values(rowOffset - 4)  = -1.0;
      values(rowOffset - 3)  = -1.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const zFaceFETag&, const ordinal_type idx) const {
    /*******************/
    /*   z == 0 face   */
    /*******************/
    // Compute row index
    ordinal_type j      = idx / (nx - 2);
    ordinal_type i      = idx % (nx - 2);
    ordinal_type rowIdx = (j + 1) * nx + i + 1;

    // Compute rowOffset
    size_type rowOffset = size_type(j) * numEntriesFrontRow + numEntriesBottomFrontRow +
                          size_type(i + 1) * faceStencilLength + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 18) = rowIdx - nx - 1;
    columns(rowOffset - 17) = rowIdx - nx;
    columns(rowOffset - 16) = rowIdx - nx + 1;
    columns(rowOffset - 15) = rowIdx - 1;
    columns(rowOffset - 14) = rowIdx;
    columns(rowOffset - 13) = rowIdx + 1;
    columns(rowOffset - 12) = rowIdx + nx - 1;
    columns(rowOffset - 11) = rowIdx + nx;
    columns(rowOffset - 10) = rowIdx + nx + 1;
    columns(rowOffset - 9)  = rowIdx + nx * ny - nx - 1;
    columns(rowOffset - 8)  = rowIdx + nx * ny - nx;
    columns(rowOffset - 7)  = rowIdx + nx * ny - nx + 1;
    columns(rowOffset - 6)  = rowIdx + nx * ny - 1;
    columns(rowOffset - 5)  = rowIdx + nx * ny;
    columns(rowOffset - 4)  = rowIdx + nx * ny + 1;
    columns(rowOffset - 3)  = rowIdx + nx * ny + nx - 1;
    columns(rowOffset - 2)  = rowIdx + nx * ny + nx;
    columns(rowOffset - 1)  = rowIdx + nx * ny + nx + 1;
    if (bottomBC == 1) {
      // Fill values
      values(rowOffset - 18) = 0.0;
      values(rowOffset - 17) = 0.0;
      values(rowOffset - 16) = 0.0;
      values(rowOffset - 15) = 0.0;
      values(rowOffset - 14) = 1.0;
      values(rowOffset - 13) = 0.0;
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 18) = -1.0;
      values(rowOffset - 17) = 0.0;
      values(rowOffset - 16) = -1.0;
      values(rowOffset - 15) = 0.0;
      values(rowOffset - 14) = 16.0;
      values(rowOffset - 13) = 0.0;
      values(rowOffset - 12) = -1.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = -1.0;
      values(rowOffset - 9)  = -1.0;
      values(rowOffset - 8)  = -2.0;
      values(rowOffset - 7)  = -1.0;
      values(rowOffset - 6)  = -2.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = -2.0;
      values(rowOffset - 3)  = -1.0;
      values(rowOffset - 2)  = -2.0;
      values(rowOffset - 1)  = -1.0;
    }

    /********************/
    /*   z == nz face   */
    /********************/
    // Compute row index
    ordinal_type k = nz - 2;
    j              = idx / (nx - 2);
    i              = idx % (nx - 2);
    rowIdx         = (k + 1) * ny * nx + (j + 1) * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesFrontRow +
                numEntriesBottomFrontRow + size_type(i + 1) * faceStencilLength + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 18) = rowIdx - nx * ny - nx - 1;
    columns(rowOffset - 17) = rowIdx - nx * ny - nx;
    columns(rowOffset - 16) = rowIdx - nx * ny - nx + 1;
    columns(rowOffset - 15) = rowIdx - nx * ny - 1;
    columns(rowOffset - 14) = rowIdx - nx * ny;
    columns(rowOffset - 13) = rowIdx - nx * ny + 1;
    columns(rowOffset - 12) = rowIdx - nx * ny + nx - 1;
    columns(rowOffset - 11) = rowIdx - nx * ny + nx;
    columns(rowOffset - 10) = rowIdx - nx * ny + nx + 1;
    columns(rowOffset - 9)  = rowIdx - nx - 1;
    columns(rowOffset - 8)  = rowIdx - nx;
    columns(rowOffset - 7)  = rowIdx - nx + 1;
    columns(rowOffset - 6)  = rowIdx - 1;
    columns(rowOffset - 5)  = rowIdx;
    columns(rowOffset - 4)  = rowIdx + 1;
    columns(rowOffset - 3)  = rowIdx + nx - 1;
    columns(rowOffset - 2)  = rowIdx + nx;
    columns(rowOffset - 1)  = rowIdx + nx + 1;
    if (topBC == 1) {
      // Fill values
      values(rowOffset - 18) = 0.0;
      values(rowOffset - 17) = 0.0;
      values(rowOffset - 16) = 0.0;
      values(rowOffset - 15) = 0.0;
      values(rowOffset - 14) = 0.0;
      values(rowOffset - 13) = 0.0;
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 1.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 18) = -1.0;
      values(rowOffset - 17) = -2.0;
      values(rowOffset - 16) = -1.0;
      values(rowOffset - 15) = -2.0;
      values(rowOffset - 14) = 0.0;
      values(rowOffset - 13) = -2.0;
      values(rowOffset - 12) = -1.0;
      values(rowOffset - 11) = -2.0;
      values(rowOffset - 10) = -1.0;
      values(rowOffset - 9)  = -1.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = -1.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 16.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = -1.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = -1.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const xEdgeFETag&, const ordinal_type idx) const {
    // Compute row index
    ordinal_type i      = idx;
    ordinal_type rowIdx = i + 1;

    // Compute rowOffset
    size_type rowOffset = size_type(i + 1) * edgeStencilLength + cornerStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - 1;
    columns(rowOffset - 11) = rowIdx;
    columns(rowOffset - 10) = rowIdx + 1;
    columns(rowOffset - 9)  = rowIdx + nx - 1;
    columns(rowOffset - 8)  = rowIdx + nx;
    columns(rowOffset - 7)  = rowIdx + nx + 1;
    columns(rowOffset - 6)  = rowIdx + nx * ny - 1;
    columns(rowOffset - 5)  = rowIdx + nx * ny;
    columns(rowOffset - 4)  = rowIdx + nx * ny + 1;
    columns(rowOffset - 3)  = rowIdx + nx * ny + nx - 1;
    columns(rowOffset - 2)  = rowIdx + nx * ny + nx;
    columns(rowOffset - 1)  = rowIdx + nx * ny + nx + 1;
    if (bottomBC == 1 || frontBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 1.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 8.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = -1.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = -1.0;
      values(rowOffset - 6)  = -1.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = -1.0;
      values(rowOffset - 3)  = -1.0;
      values(rowOffset - 2)  = -2.0;
      values(rowOffset - 1)  = -1.0;
    }

    // Compute row index
    ordinal_type j = ny - 2;
    i              = idx;
    rowIdx         = (j + 1) * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(j) * numEntriesFrontRow + numEntriesBottomFrontRow + size_type(i + 1) * edgeStencilLength +
                cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - nx - 1;
    columns(rowOffset - 11) = rowIdx - nx;
    columns(rowOffset - 10) = rowIdx - nx + 1;
    columns(rowOffset - 9)  = rowIdx - 1;
    columns(rowOffset - 8)  = rowIdx;
    columns(rowOffset - 7)  = rowIdx + 1;
    columns(rowOffset - 6)  = rowIdx + nx * ny - nx - 1;
    columns(rowOffset - 5)  = rowIdx + nx * ny - nx;
    columns(rowOffset - 4)  = rowIdx + nx * ny - nx + 1;
    columns(rowOffset - 3)  = rowIdx + nx * ny - 1;
    columns(rowOffset - 2)  = rowIdx + nx * ny;
    columns(rowOffset - 1)  = rowIdx + nx * ny + 1;
    if (bottomBC == 1 || backBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 1.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = -1.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = -1.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 8.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = -1.0;
      values(rowOffset - 5)  = -2.0;
      values(rowOffset - 4)  = -1.0;
      values(rowOffset - 3)  = -1.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = -1.0;
    }

    // Compute row index
    ordinal_type k = nz - 2;
    i              = idx;
    rowIdx         = (k + 1) * ny * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(i + 1) * edgeStencilLength +
                cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - nx * ny - 1;
    columns(rowOffset - 11) = rowIdx - nx * ny;
    columns(rowOffset - 10) = rowIdx - nx * ny + 1;
    columns(rowOffset - 9)  = rowIdx - nx * ny + nx - 1;
    columns(rowOffset - 8)  = rowIdx - nx * ny + nx;
    columns(rowOffset - 7)  = rowIdx - nx * ny + nx + 1;
    columns(rowOffset - 6)  = rowIdx - 1;
    columns(rowOffset - 5)  = rowIdx;
    columns(rowOffset - 4)  = rowIdx + 1;
    columns(rowOffset - 3)  = rowIdx + nx - 1;
    columns(rowOffset - 2)  = rowIdx + nx;
    columns(rowOffset - 1)  = rowIdx + nx + 1;
    if (topBC == 1 || frontBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 1.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = -1.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = -1.0;
      values(rowOffset - 9)  = -1.0;
      values(rowOffset - 8)  = -2.0;
      values(rowOffset - 7)  = -1.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 8.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = -1.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = -1.0;
    }

    // Compute row index
    k      = nz - 2;
    j      = ny - 2;
    i      = idx;
    rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesFrontRow +
                numEntriesBottomFrontRow + size_type(i + 1) * edgeStencilLength + cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - nx * ny - nx - 1;
    columns(rowOffset - 11) = rowIdx - nx * ny - nx;
    columns(rowOffset - 10) = rowIdx - nx * ny - nx + 1;
    columns(rowOffset - 9)  = rowIdx - nx * ny - 1;
    columns(rowOffset - 8)  = rowIdx - nx * ny;
    columns(rowOffset - 7)  = rowIdx - nx * ny + 1;
    columns(rowOffset - 6)  = rowIdx - nx - 1;
    columns(rowOffset - 5)  = rowIdx - nx;
    columns(rowOffset - 4)  = rowIdx - nx + 1;
    columns(rowOffset - 3)  = rowIdx - 1;
    columns(rowOffset - 2)  = rowIdx;
    columns(rowOffset - 1)  = rowIdx + 1;
    if (topBC == 1 || backBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 1.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = -1.0;
      values(rowOffset - 11) = -2.0;
      values(rowOffset - 10) = -1.0;
      values(rowOffset - 9)  = -1.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = -1.0;
      values(rowOffset - 6)  = -1.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = -1.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 8.0;
      values(rowOffset - 1)  = 0.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const yEdgeFETag&, const ordinal_type idx) const {
    // Compute row index
    ordinal_type j      = idx;
    ordinal_type rowIdx = (j + 1) * nx;

    // Compute rowOffset
    size_type rowOffset = size_type(j) * numEntriesFrontRow + numEntriesBottomFrontRow + edgeStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - nx;
    columns(rowOffset - 11) = rowIdx - nx + 1;
    columns(rowOffset - 10) = rowIdx;
    columns(rowOffset - 9)  = rowIdx + 1;
    columns(rowOffset - 8)  = rowIdx + nx;
    columns(rowOffset - 7)  = rowIdx + nx + 1;
    columns(rowOffset - 6)  = rowIdx + nx * ny - nx;
    columns(rowOffset - 5)  = rowIdx + nx * ny - nx + 1;
    columns(rowOffset - 4)  = rowIdx + nx * ny;
    columns(rowOffset - 3)  = rowIdx + nx * ny + 1;
    columns(rowOffset - 2)  = rowIdx + nx * ny + nx;
    columns(rowOffset - 1)  = rowIdx + nx * ny + nx + 1;
    if (bottomBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 1.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = -1.0;
      values(rowOffset - 10) = 8.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = -1.0;
      values(rowOffset - 6)  = -1.0;
      values(rowOffset - 5)  = -1.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = -2.0;
      values(rowOffset - 2)  = -1.0;
      values(rowOffset - 1)  = -1.0;
    }

    // Compute row index
    j              = idx;
    ordinal_type i = nx - 1;
    rowIdx         = (j + 1) * nx + i;

    // Compute rowOffset
    rowOffset          = size_type(j + 1) * numEntriesFrontRow + numEntriesBottomFrontRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - nx - 1;
    columns(rowOffset - 11) = rowIdx - nx;
    columns(rowOffset - 10) = rowIdx - 1;
    columns(rowOffset - 9)  = rowIdx;
    columns(rowOffset - 8)  = rowIdx + nx - 1;
    columns(rowOffset - 7)  = rowIdx + nx;
    columns(rowOffset - 6)  = rowIdx + nx * ny - nx - 1;
    columns(rowOffset - 5)  = rowIdx + nx * ny - nx;
    columns(rowOffset - 4)  = rowIdx + nx * ny - 1;
    columns(rowOffset - 3)  = rowIdx + nx * ny;
    columns(rowOffset - 2)  = rowIdx + nx * ny + nx - 1;
    columns(rowOffset - 1)  = rowIdx + nx * ny + nx;
    if (bottomBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 1.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = -1.0;
      values(rowOffset - 11) = -1.0;
      values(rowOffset - 10) = -2.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = -1.0;
      values(rowOffset - 7)  = -1.0;
      values(rowOffset - 6)  = -1.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 8.0;
      values(rowOffset - 2)  = -1.0;
      values(rowOffset - 1)  = 0.0;
    }

    // Compute row index
    ordinal_type k = nz - 2;
    j              = idx;
    rowIdx         = (k + 1) * ny * nx + (j + 1) * nx;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesFrontRow +
                numEntriesBottomFrontRow + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - nx * ny - 1;
    columns(rowOffset - 11) = rowIdx - nx * ny;
    columns(rowOffset - 10) = rowIdx - nx * ny + 1;
    columns(rowOffset - 9)  = rowIdx - nx * ny + nx - 1;
    columns(rowOffset - 8)  = rowIdx - nx * ny + nx;
    columns(rowOffset - 7)  = rowIdx - nx * ny + nx + 1;
    columns(rowOffset - 6)  = rowIdx - 1;
    columns(rowOffset - 5)  = rowIdx;
    columns(rowOffset - 4)  = rowIdx + 1;
    columns(rowOffset - 3)  = rowIdx + nx - 1;
    columns(rowOffset - 2)  = rowIdx + nx;
    columns(rowOffset - 1)  = rowIdx + nx + 1;
    if (topBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 1.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = -1.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = -1.0;
      values(rowOffset - 9)  = -1.0;
      values(rowOffset - 8)  = -2.0;
      values(rowOffset - 7)  = -1.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 8.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = -1.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = -1.0;
    }

    // Compute row index
    k      = nz - 2;
    j      = idx;
    i      = nx - 1;
    rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j + 1) * numEntriesFrontRow +
                numEntriesBottomFrontRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - nx * ny - nx - 1;
    columns(rowOffset - 11) = rowIdx - nx * ny - nx;
    columns(rowOffset - 10) = rowIdx - nx * ny - 1;
    columns(rowOffset - 9)  = rowIdx - nx * ny;
    columns(rowOffset - 8)  = rowIdx - nx * ny + nx - 1;
    columns(rowOffset - 7)  = rowIdx - nx * ny + nx;
    columns(rowOffset - 6)  = rowIdx - nx - 1;
    columns(rowOffset - 5)  = rowIdx - nx;
    columns(rowOffset - 4)  = rowIdx - 1;
    columns(rowOffset - 3)  = rowIdx;
    columns(rowOffset - 2)  = rowIdx + nx - 1;
    columns(rowOffset - 1)  = rowIdx + nx;
    if (topBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 1.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = -1.0;
      values(rowOffset - 11) = -1.0;
      values(rowOffset - 10) = -2.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = -1.0;
      values(rowOffset - 7)  = -1.0;
      values(rowOffset - 6)  = -1.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 8.0;
      values(rowOffset - 2)  = -1.0;
      values(rowOffset - 1)  = 0.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const zEdgeFETag&, const ordinal_type idx) const {
    // Compute row index
    ordinal_type k      = idx;
    ordinal_type rowIdx = (k + 1) * ny * nx;

    // Compute rowOffset
    size_type rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + edgeStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - nx * ny;
    columns(rowOffset - 11) = rowIdx - nx * ny + 1;
    columns(rowOffset - 10) = rowIdx - nx * ny + nx;
    columns(rowOffset - 9)  = rowIdx - nx * ny + nx + 1;
    columns(rowOffset - 8)  = rowIdx;
    columns(rowOffset - 7)  = rowIdx + 1;
    columns(rowOffset - 6)  = rowIdx + nx;
    columns(rowOffset - 5)  = rowIdx + nx + 1;
    columns(rowOffset - 4)  = rowIdx + ny * nx;
    columns(rowOffset - 3)  = rowIdx + ny * nx + 1;
    columns(rowOffset - 2)  = rowIdx + ny * nx + nx;
    columns(rowOffset - 1)  = rowIdx + ny * nx + nx + 1;
    if (frontBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 1.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = -1.0;
      values(rowOffset - 10) = -1.0;
      values(rowOffset - 9)  = -1.0;
      values(rowOffset - 8)  = 8.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = -2.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = -1.0;
      values(rowOffset - 2)  = -1.0;
      values(rowOffset - 1)  = -1.0;
    }

    // Compute row index
    k              = idx;
    ordinal_type i = nx - 2;
    rowIdx         = (k + 1) * ny * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(i) * faceStencilLength +
                2 * edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - nx * ny - 1;
    columns(rowOffset - 11) = rowIdx - nx * ny;
    columns(rowOffset - 10) = rowIdx - nx * ny + nx - 1;
    columns(rowOffset - 9)  = rowIdx - nx * ny + nx;
    columns(rowOffset - 8)  = rowIdx - 1;
    columns(rowOffset - 7)  = rowIdx;
    columns(rowOffset - 6)  = rowIdx + nx - 1;
    columns(rowOffset - 5)  = rowIdx + nx;
    columns(rowOffset - 4)  = rowIdx + ny * nx - 1;
    columns(rowOffset - 3)  = rowIdx + ny * nx;
    columns(rowOffset - 2)  = rowIdx + ny * nx + nx - 1;
    columns(rowOffset - 1)  = rowIdx + ny * nx + nx;
    if (frontBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 1.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = -1.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = -1.0;
      values(rowOffset - 9)  = -1.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 8.0;
      values(rowOffset - 6)  = -2.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = -1.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = -1.0;
      values(rowOffset - 1)  = -1.0;
    }

    // Compute row index
    k              = idx;
    ordinal_type j = ny - 2;
    rowIdx         = (k + 1) * ny * nx + (j + 1) * nx;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesPerGridRow +
                numEntriesFrontRow + edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - nx * ny - nx;
    columns(rowOffset - 11) = rowIdx - nx * ny - nx + 1;
    columns(rowOffset - 10) = rowIdx - nx * ny;
    columns(rowOffset - 9)  = rowIdx - nx * ny + 1;
    columns(rowOffset - 8)  = rowIdx - nx;
    columns(rowOffset - 7)  = rowIdx - nx + 1;
    columns(rowOffset - 6)  = rowIdx;
    columns(rowOffset - 5)  = rowIdx + 1;
    columns(rowOffset - 4)  = rowIdx + ny * nx - nx;
    columns(rowOffset - 3)  = rowIdx + ny * nx - nx + 1;
    columns(rowOffset - 2)  = rowIdx + ny * nx;
    columns(rowOffset - 1)  = rowIdx + ny * nx + 1;
    if (backBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 1.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = -1.0;
      values(rowOffset - 11) = -1.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = -1.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = -2.0;
      values(rowOffset - 6)  = 8.0;
      values(rowOffset - 5)  = 0.0;
      values(rowOffset - 4)  = -1.0;
      values(rowOffset - 3)  = -1.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = -1.0;
    }

    // Compute row index
    k      = idx;
    j      = ny - 2;
    i      = nx - 2;
    rowIdx = (k + 1) * ny * nx + (j + 1) * nx + i + 1;

    // Compute rowOffset
    rowOffset = size_type(k) * numEntriesPerGridPlane + numEntriesBottomPlane + size_type(j) * numEntriesPerGridRow +
                numEntriesFrontRow + size_type(i) * faceStencilLength + 2 * edgeStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 12) = rowIdx - nx * ny - nx - 1;
    columns(rowOffset - 11) = rowIdx - nx * ny - nx;
    columns(rowOffset - 10) = rowIdx - nx * ny - 1;
    columns(rowOffset - 9)  = rowIdx - nx * ny;
    columns(rowOffset - 8)  = rowIdx - nx - 1;
    columns(rowOffset - 7)  = rowIdx - nx;
    columns(rowOffset - 6)  = rowIdx - 1;
    columns(rowOffset - 5)  = rowIdx;
    columns(rowOffset - 4)  = rowIdx + ny * nx - nx - 1;
    columns(rowOffset - 3)  = rowIdx + ny * nx - nx;
    columns(rowOffset - 2)  = rowIdx + ny * nx - 1;
    columns(rowOffset - 1)  = rowIdx + ny * nx;
    if (backBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 12) = 0.0;
      values(rowOffset - 11) = 0.0;
      values(rowOffset - 10) = 0.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = 0.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 1.0;
      values(rowOffset - 4)  = 0.0;
      values(rowOffset - 3)  = 0.0;
      values(rowOffset - 2)  = 0.0;
      values(rowOffset - 1)  = 0.0;
    } else {
      // Fill values
      values(rowOffset - 12) = -1.0;
      values(rowOffset - 11) = -1.0;
      values(rowOffset - 10) = -1.0;
      values(rowOffset - 9)  = 0.0;
      values(rowOffset - 8)  = -2.0;
      values(rowOffset - 7)  = 0.0;
      values(rowOffset - 6)  = 0.0;
      values(rowOffset - 5)  = 8.0;
      values(rowOffset - 4)  = -1.0;
      values(rowOffset - 3)  = -1.0;
      values(rowOffset - 2)  = -1.0;
      values(rowOffset - 1)  = 0.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const cornerFETag&, const ordinal_type /*idx*/) const {
    // Bottom corners
    ordinal_type rowIdx = 0;
    size_type rowOffset = cornerStencilLength;
    rowmap(rowIdx + 1)  = rowOffset;

    // Fill column indices
    columns(rowOffset - 8) = rowIdx;
    columns(rowOffset - 7) = rowIdx + 1;
    columns(rowOffset - 6) = rowIdx + nx;
    columns(rowOffset - 5) = rowIdx + nx + 1;
    columns(rowOffset - 4) = rowIdx + ny * nx;
    columns(rowOffset - 3) = rowIdx + ny * nx + 1;
    columns(rowOffset - 2) = rowIdx + ny * nx + nx;
    columns(rowOffset - 1) = rowIdx + ny * nx + nx + 1;
    if (bottomBC == 1 || frontBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 8) = 1.0;
      values(rowOffset - 7) = 0.0;
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 8) = 4.0;
      values(rowOffset - 7) = 0.0;
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    rowIdx             = nx - 1;
    rowOffset          = numEntriesBottomFrontRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 8) = rowIdx - 1;
    columns(rowOffset - 7) = rowIdx;
    columns(rowOffset - 6) = rowIdx + nx - 1;
    columns(rowOffset - 5) = rowIdx + nx;
    columns(rowOffset - 4) = rowIdx + ny * nx - 1;
    columns(rowOffset - 3) = rowIdx + ny * nx;
    columns(rowOffset - 2) = rowIdx + ny * nx + nx - 1;
    columns(rowOffset - 1) = rowIdx + ny * nx + nx;
    if (bottomBC == 1 || frontBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 8) = 0.0;
      values(rowOffset - 7) = 1.0;
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 8) = 0.0;
      values(rowOffset - 7) = 4.0;
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = -1.0;
    }

    rowIdx             = (ny - 1) * nx;
    rowOffset          = size_type(ny - 2) * numEntriesFrontRow + numEntriesBottomFrontRow + cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 8) = rowIdx - nx;
    columns(rowOffset - 7) = rowIdx - nx + 1;
    columns(rowOffset - 6) = rowIdx;
    columns(rowOffset - 5) = rowIdx + 1;
    columns(rowOffset - 4) = rowIdx + ny * nx - nx;
    columns(rowOffset - 3) = rowIdx + ny * nx - nx + 1;
    columns(rowOffset - 2) = rowIdx + ny * nx;
    columns(rowOffset - 1) = rowIdx + ny * nx + 1;
    if (bottomBC == 1 || backBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 8) = 0.0;
      values(rowOffset - 7) = 0.0;
      values(rowOffset - 6) = 1.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 8) = 0.0;
      values(rowOffset - 7) = -1.0;
      values(rowOffset - 6) = 4.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = -1.0;
    }

    rowIdx             = ny * nx - 1;
    rowOffset          = numEntriesBottomPlane;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 8) = rowIdx - nx - 1;
    columns(rowOffset - 7) = rowIdx - nx;
    columns(rowOffset - 6) = rowIdx - 1;
    columns(rowOffset - 5) = rowIdx;
    columns(rowOffset - 4) = rowIdx + ny * nx - nx - 1;
    columns(rowOffset - 3) = rowIdx + ny * nx - nx;
    columns(rowOffset - 2) = rowIdx + ny * nx - 1;
    columns(rowOffset - 1) = rowIdx + ny * nx;
    if (bottomBC == 1 || backBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 8) = 0.0;
      values(rowOffset - 7) = 0.0;
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 1.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 8) = -1.0;
      values(rowOffset - 7) = 0.0;
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 4.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = 0.0;
    }

    // Top corners
    rowIdx             = (nz - 1) * ny * nx;
    rowOffset          = size_type(nz - 2) * numEntriesPerGridPlane + numEntriesBottomPlane + cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 8) = rowIdx - ny * nx;
    columns(rowOffset - 7) = rowIdx - ny * nx + 1;
    columns(rowOffset - 6) = rowIdx - ny * nx + nx;
    columns(rowOffset - 5) = rowIdx - ny * nx + nx + 1;
    columns(rowOffset - 4) = rowIdx;
    columns(rowOffset - 3) = rowIdx + 1;
    columns(rowOffset - 2) = rowIdx + nx;
    columns(rowOffset - 1) = rowIdx + nx + 1;
    if (topBC == 1 || frontBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 8) = 0.0;
      values(rowOffset - 7) = 0.0;
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 8) = 0.0;
      values(rowOffset - 7) = -1.0;
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = 4.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = -1.0;
    }

    rowIdx             = (nz - 1) * ny * nx + nx - 1;
    rowOffset          = size_type(nz - 2) * numEntriesPerGridPlane + numEntriesBottomPlane + numEntriesBottomFrontRow;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 8) = rowIdx - ny * nx - 1;
    columns(rowOffset - 7) = rowIdx - ny * nx;
    columns(rowOffset - 6) = rowIdx - ny * nx + nx - 1;
    columns(rowOffset - 5) = rowIdx - ny * nx + nx;
    columns(rowOffset - 4) = rowIdx - 1;
    columns(rowOffset - 3) = rowIdx;
    columns(rowOffset - 2) = rowIdx + nx - 1;
    columns(rowOffset - 1) = rowIdx + nx;
    if (topBC == 1 || frontBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 8) = 0.0;
      values(rowOffset - 7) = 0.0;
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 1.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 8) = -1.0;
      values(rowOffset - 7) = 0.0;
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 4.0;
      values(rowOffset - 2) = -1.0;
      values(rowOffset - 1) = 0.0;
    }

    rowIdx    = nz * ny * nx - nx;
    rowOffset = size_type(nz - 2) * numEntriesPerGridPlane + numEntriesBottomPlane +
                size_type(ny - 2) * numEntriesFrontRow + numEntriesBottomFrontRow + cornerStencilLength;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 8) = rowIdx - ny * nx - nx;
    columns(rowOffset - 7) = rowIdx - ny * nx - nx + 1;
    columns(rowOffset - 6) = rowIdx - ny * nx;
    columns(rowOffset - 5) = rowIdx - ny * nx + 1;
    columns(rowOffset - 4) = rowIdx - nx;
    columns(rowOffset - 3) = rowIdx - nx + 1;
    columns(rowOffset - 2) = rowIdx;
    columns(rowOffset - 1) = rowIdx + 1;
    if (topBC == 1 || backBC == 1 || leftBC == 1) {
      // Fill values
      values(rowOffset - 8) = 0.0;
      values(rowOffset - 7) = 0.0;
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 1.0;
      values(rowOffset - 1) = 0.0;
    } else {
      // Fill values
      values(rowOffset - 8) = -1.0;
      values(rowOffset - 7) = -1.0;
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = -1.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = -1.0;
      values(rowOffset - 2) = 4.0;
      values(rowOffset - 1) = 0.0;
    }

    rowIdx             = nz * ny * nx - 1;
    rowOffset          = numEntries;
    rowmap(rowIdx + 1) = rowOffset;

    // Fill column indices
    columns(rowOffset - 8) = rowIdx - ny * nx - nx - 1;
    columns(rowOffset - 7) = rowIdx - ny * nx - nx;
    columns(rowOffset - 6) = rowIdx - ny * nx - 1;
    columns(rowOffset - 5) = rowIdx - ny * nx;
    columns(rowOffset - 4) = rowIdx - nx - 1;
    columns(rowOffset - 3) = rowIdx - nx;
    columns(rowOffset - 2) = rowIdx - 1;
    columns(rowOffset - 1) = rowIdx;
    if (topBC == 1 || backBC == 1 || rightBC == 1) {
      // Fill values
      values(rowOffset - 8) = 0.0;
      values(rowOffset - 7) = 0.0;
      values(rowOffset - 6) = 0.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = 0.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 1.0;
    } else {
      // Fill values
      values(rowOffset - 8) = -1.0;
      values(rowOffset - 7) = -1.0;
      values(rowOffset - 6) = -1.0;
      values(rowOffset - 5) = 0.0;
      values(rowOffset - 4) = -1.0;
      values(rowOffset - 3) = 0.0;
      values(rowOffset - 2) = 0.0;
      values(rowOffset - 1) = 4.0;
    }
  }
};

template <typename CrsMatrix_t, typename mat_structure>
CrsMatrix_t generate_structured_matrix3D(const std::string stencil, const mat_structure& structure) {
  typedef typename CrsMatrix_t::StaticCrsGraphType graph_t;
  typedef typename CrsMatrix_t::row_map_type::non_const_type row_map_view_t;
  typedef typename CrsMatrix_t::index_type::non_const_type cols_view_t;
  typedef typename CrsMatrix_t::values_type::non_const_type scalar_view_t;
  typedef typename CrsMatrix_t::non_const_size_type size_type;
  typedef typename CrsMatrix_t::non_const_ordinal_type ordinal_type;

  int stencil_type = 0;
  if (stencil == "FD") {
    stencil_type = FD;
  } else if (stencil == "FE") {
    stencil_type = FE;
  } else {
    std::ostringstream os;
    os << "Test::generate_structured_matrix3D only accepts stencil: FD and "
          "FEM, you passed: "
       << stencil << " !" << std::endl;
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Extract geometric data
  const ordinal_type nx          = structure(0, 0);
  const ordinal_type ny          = structure(1, 0);
  const ordinal_type nz          = structure(2, 0);
  const ordinal_type numNodes    = ny * nx * nz;
  const ordinal_type leftBC      = structure(0, 1);
  const ordinal_type rightBC     = structure(0, 2);
  const ordinal_type frontBC     = structure(1, 1);
  const ordinal_type backBC      = structure(1, 2);
  const ordinal_type bottomBC    = structure(2, 1);
  const ordinal_type topBC       = structure(2, 2);
  const ordinal_type numInterior = (nx - leftBC - rightBC) * (ny - frontBC - backBC) * (nz - bottomBC - topBC);
  const ordinal_type numFace     = (leftBC + rightBC) * (ny - frontBC - backBC) * (nz - bottomBC - topBC) +
                               (frontBC + backBC) * (nx - leftBC - rightBC) * (nz - bottomBC - topBC) +
                               (bottomBC + topBC) * (nx - leftBC - rightBC) * (ny - frontBC - backBC);
  const ordinal_type numEdge =
      (frontBC * bottomBC + frontBC * topBC + backBC * bottomBC + backBC * topBC) * (nx - leftBC - rightBC) +
      (leftBC * bottomBC + leftBC * topBC + rightBC * bottomBC + rightBC * topBC) * (ny - frontBC - backBC) +
      (leftBC * frontBC + leftBC * backBC + rightBC * frontBC + rightBC * backBC) * (nz - bottomBC - topBC);
  const ordinal_type numCorner = leftBC * frontBC * bottomBC + rightBC * frontBC * bottomBC +
                                 leftBC * backBC * bottomBC + rightBC * backBC * bottomBC + leftBC * frontBC * topBC +
                                 rightBC * frontBC * topBC + leftBC * backBC * topBC + rightBC * backBC * topBC;
  ordinal_type interiorStencilLength = 0, faceStencilLength = 0, edgeStencilLength = 0, cornerStencilLength = 0;

  if (stencil_type == FD) {
    interiorStencilLength = 7;
    faceStencilLength     = 6;
    edgeStencilLength     = 5;
    cornerStencilLength   = 4;
  } else if (stencil_type == FE) {
    interiorStencilLength = 27;
    faceStencilLength     = 18;
    edgeStencilLength     = 12;
    cornerStencilLength   = 8;
  }

  const size_type numEntries = numInterior * interiorStencilLength + numFace * faceStencilLength +
                               numEdge * edgeStencilLength + numCorner * cornerStencilLength;

  // Create matrix data
  row_map_view_t rowmap_view("rowmap_view", numNodes + 1);
  cols_view_t columns_view("colsmap_view", numEntries);
  scalar_view_t values_view("values_view", numEntries);

  // Fill the CrsGraph and the CrsMatrix
  // To start simple we construct 2D 5pt stencil Laplacian.
  // We assume Neumann boundary conditions on the edge of the domain.

  fill_3D_matrix_functor<CrsMatrix_t> fill_3D_matrix(stencil_type, nx, ny, nz, leftBC, rightBC, frontBC, backBC,
                                                     bottomBC, topBC, rowmap_view, columns_view, values_view);

  fill_3D_matrix.compute();

  graph_t static_graph(columns_view, rowmap_view);
  std::string name;
  if (stencil_type == FD) {
    name = "CrsMatrixFD";
  } else if (stencil_type == FE) {
    name = "CrsMatrixFE";
  }

  return CrsMatrix_t(name, numNodes, values_view, static_graph);

}  // generate_structured_matrix3D

}  // namespace Test

#endif  // KOKKOSKERNELS_TEST_STRUCTURE_MATRIX_HPP
