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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_CubatureControlVolumeSideDef.hpp
    \brief  Header file for the Intrepid2::CubatureControlVolume class.
    \author Created by K. Peterson, P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_CONTROLVOLUME_SIDE_DEF_HPP__
#define __INTREPID2_CUBATURE_CONTROLVOLUME_SIDE_DEF_HPP__

namespace Intrepid2 {

  template <typename SpT, typename PT, typename WT>
  CubatureControlVolumeSide<SpT,PT,WT>::
  CubatureControlVolumeSide(const shards::CellTopology cellTopology) {

    // define primary cell topology with given one
    primaryCellTopo_ = cellTopology;

    // subcell is defined either hex or quad according to dimension
    const auto spaceDim = primaryCellTopo_.getDimension();
    switch (spaceDim) {
    case 2:
      subcvCellTopo_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >());
      break;
    case 3:
      subcvCellTopo_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >());
      break;
    }

    // computation order is always one;
    degree_ = 1;

    // precompute reference side points that are repeatedly used in get cubature
    const auto maxNumNodesPerSide = 10;
    const auto numSideNodeMaps = (spaceDim == 2 ? 1 : 2);

    const ordinal_type sideOrd[2] = { 1, 5 };
    Kokkos::View<ordinal_type**,Kokkos::HostSpace> sideNodeMapHost("CubatureControlVolumeSide::sideNodeMapHost",
                                                                   numSideNodeMaps, maxNumNodesPerSide);

    const auto sideDim = spaceDim - 1;
    for (ordinal_type i=0;i<numSideNodeMaps;++i) {
      const auto side = sideOrd[i];
      sideNodeMapHost(i,0) = subcvCellTopo_.getNodeCount(sideDim, side);
      for (ordinal_type j=0;j<sideNodeMapHost(i,0);++j)
        sideNodeMapHost(i,j+1) = subcvCellTopo_.getNodeMap(sideDim, side, j);
    }
    sideNodeMap_ = Kokkos::create_mirror_view(typename SpT::memory_space(), sideNodeMapHost);
    Kokkos::deep_copy(sideNodeMap_, sideNodeMapHost);

    Kokkos::DynRankView<PT,SpT> sideCenterLocal("CubatureControlVolumeSide::sideCenterLocal",
                                                1, sideDim);

    // map to reference subcell function relies on uvm; some utility functions in cell tools still need uvm
    sidePoints_ = Kokkos::DynRankView<PT,SpT>("CubatureControlVolumeSide::sidePoints", numSideNodeMaps, spaceDim);
    for (ordinal_type i=0;i<numSideNodeMaps;++i) {
      const auto sideRange = Kokkos::pair<ordinal_type,ordinal_type>(i, i+1);
      auto sidePoint = Kokkos::subdynrankview(sidePoints_, sideRange, Kokkos::ALL());
      CellTools<SpT>::mapToReferenceSubcell(sidePoint,
                                            sideCenterLocal,
                                            sideDim,
                                            sideOrd[i],
                                            subcvCellTopo_);
    }
  }
  
  template <typename SpT, typename PT, typename WT>
  void
  CubatureControlVolumeSide<SpT,PT,WT>::
  getCubature( pointViewType  cubPoints,
               weightViewType cubWeights,
               pointViewType  cellCoords ) const {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolumeSide): cubPoints must have rank 3 (C,P,D).");
    INTREPID2_TEST_FOR_EXCEPTION( cubWeights.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolumeSide): cubWeights must have rank 2 (C,P,D).");
    INTREPID2_TEST_FOR_EXCEPTION( cellCoords.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolumeSide): cellCoords must have rank 3 of (C,P,D).");

    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.extent(0) != cellCoords.extent(0) ||
                                  cubPoints.extent(0) != cubWeights.extent(0), std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolumeSide): cubPoints, cubWeights and cellCoords dimension(0) are not consistent, numCells");

    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.extent(2) != cellCoords.extent(2) ||
                                  cubPoints.extent(2) != cubWeights.extent(2) ||
                                  static_cast<ordinal_type>(cubPoints.extent(2)) != getDimension(), std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolumeSide): cubPoints, cellCoords, this->getDimension() are not consistent, spaceDim.");
#endif
    typedef Kokkos::DynRankView<PT,SpT> tempPointViewType;

    // get array dimensions
    const auto numCells = cellCoords.extent(0);
    const auto numNodesPerCell = cellCoords.extent(1);
    const auto spaceDim = cellCoords.extent(2);

    const auto numNodesPerSubcv = subcvCellTopo_.getNodeCount();
    tempPointViewType subcvCoords("CubatureControlVolumeSide::subcvCoords",
                                  numCells, numNodesPerCell, numNodesPerSubcv, spaceDim);
    CellTools<SpT>::getSubcvCoords(subcvCoords,
                                   cellCoords,
                                   primaryCellTopo_);

    const auto numSideNodeMaps = (spaceDim == 2 ? 1 : 2);
    const ordinal_type sideOrd[2] = { 1, 5 };

    Kokkos::pair<ordinal_type,ordinal_type> nodeRangePerSide[2];

    // the second rage is cell specific to handle remained sides
    switch (primaryCellTopo_.getKey()) {
    case shards::Triangle<3>::key:
    case shards::Quadrilateral<4>::key:
      nodeRangePerSide[0].first  = 0;
      nodeRangePerSide[0].second = nodeRangePerSide[0].first + numNodesPerCell;
      break;
    case shards::Hexahedron<8>::key:
      nodeRangePerSide[0].first  = 0;
      nodeRangePerSide[0].second = nodeRangePerSide[0].first + numNodesPerCell;
      nodeRangePerSide[1].first  = numNodesPerCell;
      nodeRangePerSide[1].second = nodeRangePerSide[1].first + 4;
      break;
    case shards::Tetrahedron<4>::key:
      nodeRangePerSide[0].first  = 0;
      nodeRangePerSide[0].second = nodeRangePerSide[0].first + numNodesPerCell;
      nodeRangePerSide[1].first  = 3;
      nodeRangePerSide[1].second = nodeRangePerSide[1].first + 3;
      break;
    }

    for (ordinal_type i=0;i<numSideNodeMaps;++i) {
      const auto numSubcvPoints = 1;
      const auto numNodesPerThisSide = nodeRangePerSide[i].second - nodeRangePerSide[i].first;
      tempPointViewType subcvJacobian("CubatureControlVolume::subcvJacobian",
                                      numCells, numNodesPerThisSide, numSubcvPoints, spaceDim, spaceDim);
      
      tempPointViewType subcvSideNormal("CubatureControlVolume::subcvSideNormal",
                                        numCells, numNodesPerThisSide, numSubcvPoints, spaceDim);
      
      // numNodesPerCell is maximum 8; this repeated run is necessary because of cell tools input consideration
      const auto sideRange = Kokkos::pair<ordinal_type,ordinal_type>(i, i+1);
      const auto sidePoint = Kokkos::subdynrankview(sidePoints_, sideRange, Kokkos::ALL());

      for (auto node=0;node<numNodesPerThisSide;++node) {
        auto subcvJacobianNode    = Kokkos::subdynrankview(subcvJacobian,   Kokkos::ALL(), node, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto subcvCoordsNode      = Kokkos::subdynrankview(subcvCoords,     Kokkos::ALL(), node, Kokkos::ALL(), Kokkos::ALL());
        auto subcvSideNormalNode  = Kokkos::subdynrankview(subcvSideNormal, Kokkos::ALL(), node, Kokkos::ALL(), Kokkos::ALL());      

        CellTools<SpT>::setJacobian(subcvJacobianNode,    // C, P, D, D
                                    sidePoint,            //    P, D
                                    subcvCoordsNode,      // C, N, D
                                    subcvCellTopo_);
        
        CellTools<SpT>::getPhysicalSideNormals(subcvSideNormalNode,  // C, P, D
                                               subcvJacobianNode,
                                               sideOrd[i],
                                               subcvCellTopo_);    // C, P, D, D
      }
      
      typedef Kokkos::View<ordinal_type*,SpT> mapViewType;
      const auto sideNodeMap = Kokkos::subview(sideNodeMap_, i, Kokkos::ALL());
      
      typedef typename ExecSpace<typename pointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;
    
      const auto loopSize = numCells;
      Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
      
      typedef Functor<pointViewType,weightViewType,tempPointViewType,tempPointViewType,mapViewType> FunctorType;

      auto cubPointsThisSide  = Kokkos::subdynrankview(cubPoints,  Kokkos::ALL(), nodeRangePerSide[i], Kokkos::ALL());
      auto cubWeightsThisSide = Kokkos::subdynrankview(cubWeights, Kokkos::ALL(), nodeRangePerSide[i], Kokkos::ALL());

      Kokkos::parallel_for( policy, FunctorType(cubPointsThisSide,  
                                                cubWeightsThisSide, 
                                                subcvCoords,
                                                subcvSideNormal,
                                                sideNodeMap) );
    } 
  }
    
} 

#endif

