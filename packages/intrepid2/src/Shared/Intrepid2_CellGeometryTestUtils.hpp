// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_CellGeometryTestUtils.hpp
    \brief  Utility methods for working with CellGeometry objects in unit tests.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_CellGeometryTestUtils_h
#define Intrepid2_CellGeometryTestUtils_h

#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_ScalarView.hpp"

#include "Intrepid2_TestUtils.hpp"

namespace Intrepid2
{
/** \brief Use the cell nodes provided by one cell geometry object to create another CellGeometry that is node-based (as opposed to tensor-grid-based, uniform-grid-based, etc.)
   \param [in] anyCellGeometry - the projected geometry degrees of freedom
   \param [in] copyAffineness - if true, the resulting geometry will be marked as affine if the original geometry was, allowing reduce storage of Jacobians, etc.
   \return a representation of the same geometry, defined using cell-to-nodes and node-to-coordinates containers.
*/
  template< typename PointScalar, int spaceDim, typename DeviceType >
  inline
  CellGeometry<PointScalar, spaceDim, DeviceType> getNodalCellGeometry(CellGeometry<PointScalar, spaceDim, DeviceType> &anyCellGeometry,
                                                                       const bool &copyAffineness)
  {
    // copy the nodes from CellGeometry into a raw View
    const int numCells = anyCellGeometry.extent_int(0);
    const int numNodes = anyCellGeometry.extent_int(1);
    
    using PointScalarView = ScalarView<PointScalar, DeviceType >;
    using intView         = ScalarView<        int, DeviceType >;
    
    auto cellToNodes      = intView                          ("cell to nodes", numCells, numNodes);    // this will be a one-to-one mapping
    PointScalarView nodes = getView<PointScalar,DeviceType>("nodes", numCells * numNodes, spaceDim); // we store redundant copies of vertices for each cell

    using ExecutionSpace = typename DeviceType::execution_space;
    auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,numNodes});
    
    // "workset"
    Kokkos::parallel_for("copy cell nodes from CellGeometry", policy,
    KOKKOS_LAMBDA (const int &cellOrdinal, const int &nodeOrdinalInCell) {
      const int globalNodeOrdinal = cellOrdinal * numNodes + nodeOrdinalInCell;
      for (int d=0; d<spaceDim; d++)
      {
        nodes(globalNodeOrdinal,d) = anyCellGeometry(cellOrdinal,nodeOrdinalInCell,d);
      }
      cellToNodes(cellOrdinal,nodeOrdinalInCell) = globalNodeOrdinal;
    });
    
    ExecutionSpace().fence();
    
    const bool claimAffine  = copyAffineness && anyCellGeometry.affine();
    const auto nodeOrdering = anyCellGeometry.nodeOrderingForHypercubes();
    
    const auto cellTopology = anyCellGeometry.cellTopology();
    CellGeometry<PointScalar, spaceDim, DeviceType> nodalCellGeometry(cellTopology,cellToNodes,nodes,claimAffine,nodeOrdering);
    
    return nodalCellGeometry;
  }

/** \brief Create a uniform Cartesian mesh, with origin at 0, and domain extents and mesh widths that can be different in different coordinate dimensions.
   \param [in] domainExtents - array specifying the extent of the domain in each coordinate dimension.
   \param [in] gridCellCounts - array specifying the number of cells in each coordinate dimension.
   \return a uniform Cartesion mesh, with origin at 0, with the specified domain extents and grid cell counts.
*/
  template<class PointScalar, int spaceDim, typename DeviceType>
  inline CellGeometry<PointScalar,spaceDim,DeviceType> uniformCartesianMesh(const Kokkos::Array<PointScalar,spaceDim> &domainExtents,
                                                                            const Kokkos::Array<int,spaceDim> &gridCellCounts)
  {
    Kokkos::Array<PointScalar,spaceDim> origin;
    for (int d=0; d<spaceDim; d++)
    {
      origin[d] = 0.0;
    }
    
    using CellGeometry = ::Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType>;
    const auto NO_SUBDIVISION = CellGeometry::NO_SUBDIVISION;
    const auto HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS = CellGeometry::HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS;
    
    return CellGeometry(origin, domainExtents, gridCellCounts, NO_SUBDIVISION, HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS);
  }

/** \brief Create a uniform Cartesian mesh, with origin at 0, and the same extent and number of grid cells in each coordinate dimension.
   \param [in] domainExtent - the extent of the domain in each coordinate dimension.
   \param [in] gridCellCounts - the number of cells in each coordinate dimension.
   \return a uniform Cartesion mesh, with origin at 0, with the specified domain extent and grid cell counts.
*/
  template<class PointScalar, int spaceDim, typename DeviceType>
  inline CellGeometry<PointScalar,spaceDim,DeviceType> uniformCartesianMesh(const PointScalar &domainExtent, const int &meshWidth)
  {
    Kokkos::Array<PointScalar,spaceDim> origin;
    Kokkos::Array<PointScalar,spaceDim> domainExtents;
    Kokkos::Array<int,spaceDim> gridCellCounts;
    for (int d=0; d<spaceDim; d++)
    {
      origin[d] = 0.0;
      domainExtents[d] = domainExtent;
      gridCellCounts[d] = meshWidth;
    }
    
    using CellGeometry = ::Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType>;
    const auto NO_SUBDIVISION = CellGeometry::NO_SUBDIVISION;
    const auto HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS = CellGeometry::HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS;
    
    return CellGeometry(origin, domainExtents, gridCellCounts, NO_SUBDIVISION, HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS);
  }
} // namespace Intrepid2

#endif /* Intrepid2_CellGeometryTestUtils_h */
