// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   PAMatrixTests_Utils.hpp
    \brief  Meshes for PAMatrix tests.
    \author Nathan V. Roberts
*/

#ifndef PAMatrixTests_Utils_h
#define PAMatrixTests_Utils_h

#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_ScalarView.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <vector>

namespace Intrepid2
{
  /** \brief A conforming hypercube mesh with gridCellCounts cells in each dimension, with
             - graded spacing in each dimension (so that the Jacobians differ from cell to cell), and optionally a global
               shear (giving the Jacobians off-diagonal entries; the cells remain affine);
             - optionally, interior vertices perturbed by perturbation times the local spacing (making cells non-affine); and
             - global vertex numbers permuted (so that the orientations of the cells' edges and faces, which are determined by
               the global vertex numbers, are mixed, including rotated faces of hexahedra).
             Cells use the shards vertex ordering (HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS), as standard assembly expects.
  */
  template<class PointScalar, int spaceDim, typename DeviceType>
  CellGeometry<PointScalar,spaceDim,DeviceType> perturbedMixedOrientationMesh(const Kokkos::Array<int,spaceDim> &gridCellCounts,
                                                                              const double perturbation = 0.0,
                                                                              const double shear = 0.0,
                                                                              const unsigned seed = 1234)
  {
    shards::CellTopology cellTopo;
    if      (spaceDim == 1) cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >());
    else if (spaceDim == 2) cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >());
    else if (spaceDim == 3) cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >());
    INTREPID2_TEST_FOR_EXCEPTION(spaceDim > 3, std::invalid_argument, "perturbedMixedOrientationMesh() supports spaceDim <= 3");

    int numVertices = 1, numCells = 1;
    Kokkos::Array<int,spaceDim> vertexCounts;
    for (int d=0; d<spaceDim; d++)
    {
      vertexCounts[d] = gridCellCounts[d] + 1;
      numVertices *= vertexCounts[d];
      numCells    *= gridCellCounts[d];
    }
    const int verticesPerCell = 1 << spaceDim;

    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);

    // global vertex numbering: a random permutation of the lexicographic numbering
    std::vector<int> permutation(numVertices);
    std::iota(permutation.begin(), permutation.end(), 0);
    std::shuffle(permutation.begin(), permutation.end(), generator);

    using PointScalarView = ScalarView<PointScalar, DeviceType>;
    using IntView         = ScalarView<int,         DeviceType>;
    PointScalarView nodes("perturbed mesh nodes", numVertices, spaceDim);
    IntView cellToNodes("perturbed mesh cell to nodes", numCells, verticesPerCell);
    auto nodesHost       = Kokkos::create_mirror_view(nodes);
    auto cellToNodesHost = Kokkos::create_mirror_view(cellToNodes);

    // graded 1D coordinates: x_i = (i/n)^(1 + 0.3 (d+1))
    auto coordinate = [&gridCellCounts](const int d, const int i) {
      return std::pow(double(i) / gridCellCounts[d], 1.0 + 0.3 * (d+1));
    };
    for (int v=0; v<numVertices; v++)
    {
      int remainder = v;
      Kokkos::Array<double,spaceDim> x;
      for (int d=0; d<spaceDim; d++)
      {
        const int i = remainder % vertexCounts[d];
        remainder /= vertexCounts[d];
        const double h = coordinate(d, std::min(i+1, gridCellCounts[d])) - coordinate(d, std::max(i-1, 0));
        const bool interior = (i > 0) && (i < gridCellCounts[d]);
        x[d] = coordinate(d, i) + (interior ? 0.5 * perturbation * h * unit(generator) : 0.0);
      }
      // shear: x_d += shear * x_{d+1}
      for (int d=0; d<spaceDim; d++)
      {
        const double sheared = (d+1 < spaceDim) ? x[d] + shear * x[d+1] : x[d];
        nodesHost(permutation[v], d) = sheared;
      }
    }
    for (int cell=0; cell<numCells; cell++)
    {
      Kokkos::Array<int,spaceDim> cellIndex;
      int remainder = cell;
      for (int d=0; d<spaceDim; d++)
      {
        cellIndex[d] = remainder % gridCellCounts[d];
        remainder /= gridCellCounts[d];
      }
      for (int localVertex=0; localVertex<verticesPerCell; localVertex++)
      {
        // shards vertex order: counterclockwise in (x,y) -- (0,0), (1,0), (1,1), (0,1) -- then the same at z = 1
        const int xyVertex = localVertex % 4;
        const int xOffset = ((xyVertex == 1) || (xyVertex == 2)) ? 1 : 0;
        const int yOffset = (xyVertex >= 2) ? 1 : 0;
        const int zOffset = localVertex / 4;
        const int offsets[3] = {xOffset, yOffset, zOffset};
        int lexicographic = 0, stride = 1;
        for (int d=0; d<spaceDim; d++)
        {
          lexicographic += (cellIndex[d] + offsets[d]) * stride;
          stride *= vertexCounts[d];
        }
        cellToNodesHost(cell, localVertex) = permutation[lexicographic];
      }
    }
    Kokkos::deep_copy(nodes, nodesHost);
    Kokkos::deep_copy(cellToNodes, cellToNodesHost);

    const bool claimAffine = false;
    using Geometry = CellGeometry<PointScalar,spaceDim,DeviceType>;
    return Geometry(cellTopo, cellToNodes, nodes, claimAffine, Geometry::HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS);
  }
} // namespace Intrepid2

#endif /* PAMatrixTests_Utils_h */
