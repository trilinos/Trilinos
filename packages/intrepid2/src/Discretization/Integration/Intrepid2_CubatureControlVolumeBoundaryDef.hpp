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

/** \file   Intrepid2_CubatureControlVolumeBoundaryDef.hpp
    \brief  Definition file for the Intrepid2::CubatureControlVolumeBoundary class.
    \author Created by K. Peterson, P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_CONTROLVOLUME_BOUNDARY_DEF_HPP__
#define __INTREPID2_CUBATURE_CONTROLVOLUME_BOUNDARY_DEF_HPP__

namespace Intrepid2{


  template <typename SpT, typename PT, typename WT>
  CubatureControlVolumeBoundary<SpT,PT,WT>::
  CubatureControlVolumeBoundary(const shards::CellTopology cellTopology,
                                const ordinal_type sideIndex) {

    // define primary cell topology with given one
    primaryCellTopo_ = cellTopology;

    // subcell is defined either hex or quad according to dimension
    const ordinal_type spaceDim = primaryCellTopo_.getDimension();
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

    sideIndex_ = sideIndex;

    // precompute subcontrol volume side index corresponding to primary cell side
    const ordinal_type sideDim = spaceDim - 1;

    const ordinal_type numPrimarySides   = primaryCellTopo_.getSubcellCount(sideDim);
    const ordinal_type numPrimarySideNodes = primaryCellTopo_.getNodeCount(sideDim, sideIndex_);

    const ordinal_type numSubcvSideNodes = subcvCellTopo_.getNodeCount(sideDim, 0);

    boundarySidesHost_ = Kokkos::View<ordinal_type**,Kokkos::HostSpace>("CubatureControlVolumeBoundary::boundarySidesHost",
                                                                        numPrimarySides, numSubcvSideNodes);
    
    // tabulate node map on boundary side
    switch (primaryCellTopo_.getKey()) {
    case shards::Triangle<3>::key:
    case shards::Quadrilateral<4>::key: {
      for (ordinal_type i=0;i<numPrimarySides;++i) {
        boundarySidesHost_(i,0) = 0; 
        boundarySidesHost_(i,1) = 3;
      }
      break;
    }
    case shards::Hexahedron<8>::key: {
      // sides 0-3
      for (ordinal_type i=0;i<4;++i) {
        boundarySidesHost_(i,0) = 0; 
        boundarySidesHost_(i,1) = 3;
        boundarySidesHost_(i,2) = 3; 
        boundarySidesHost_(i,3) = 0;
      }

      // side 4
      boundarySidesHost_(4,0) = 4; 
      boundarySidesHost_(4,1) = 4;
      boundarySidesHost_(4,2) = 4; 
      boundarySidesHost_(4,3) = 4;

      // side 5
      boundarySidesHost_(5,0) = 5; 
      boundarySidesHost_(5,1) = 5;
      boundarySidesHost_(5,2) = 5; 
      boundarySidesHost_(5,3) = 5;
      break;
    }
    case shards::Tetrahedron<4>::key: {
      boundarySidesHost_(0,0) = 0; 
      boundarySidesHost_(0,1) = 3; 
      boundarySidesHost_(0,2) = 0;

      boundarySidesHost_(1,0) = 0; 
      boundarySidesHost_(1,1) = 3; 
      boundarySidesHost_(1,2) = 3;

      boundarySidesHost_(2,0) = 3; 
      boundarySidesHost_(2,1) = 4; 
      boundarySidesHost_(2,2) = 0;

      boundarySidesHost_(3,0) = 4; 
      boundarySidesHost_(3,1) = 4; 
      boundarySidesHost_(3,2) = 4;
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                    ">>> ERROR (CubatureControlVolumeBoundary: invalid cell topology.");
    } 
    }
    
    Kokkos::DynRankView<PT,SpT> sideCenterLocal("CubatureControlVolumeBoundary::sideCenterLocal",
                                                1, sideDim);
    // map to reference subcell function relies on uvm; some utility functions in cell tools still need uvm
    sidePoints_ = Kokkos::DynRankView<PT,SpT>("CubatureControlVolumeBoundary::sidePoints", 
                                              numPrimarySideNodes, spaceDim);

    for (ordinal_type i=0;i<numPrimarySideNodes;++i) {
      const ordinal_type sideOrd = boundarySidesHost_(sideIndex_,i);
      const auto sideRange = Kokkos::pair<ordinal_type,ordinal_type>(i, i+1);
      const auto sidePoint = Kokkos::subdynrankview(sidePoints_, sideRange, Kokkos::ALL());
      CellTools<SpT>::mapToReferenceSubcell(sidePoint,
                                            sideCenterLocal,
                                            sideDim,
                                            sideOrd,
                                            subcvCellTopo_);
    }

    const ordinal_type maxNumNodesPerSide = 10;
    Kokkos::View<ordinal_type**,Kokkos::HostSpace> sideNodeMapHost("CubatureControlVolumeSide::sideNodeMapHost",
                                                                   numPrimarySideNodes, maxNumNodesPerSide);
    for (ordinal_type i=0;i<numPrimarySideNodes;++i) {
      const ordinal_type sideOrd = boundarySidesHost_(sideIndex_,i);
      sideNodeMapHost(i,0) = subcvCellTopo_.getNodeCount(sideDim, sideOrd);
      for (ordinal_type j=0;j<sideNodeMapHost(i,0);++j)
        sideNodeMapHost(i,j+1) = subcvCellTopo_.getNodeMap(sideDim, sideOrd, j);
    }
    sideNodeMap_ = Kokkos::create_mirror_view(typename SpT::memory_space(), sideNodeMapHost);
    Kokkos::deep_copy(sideNodeMap_, sideNodeMapHost);
  }

  template <typename SpT, typename PT, typename WT>
  void
  CubatureControlVolumeBoundary<SpT,PT,WT>::
  getCubature( pointViewType  cubPoints,
               weightViewType cubWeights,
               pointViewType  cellCoords ) const {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolumeBoundary): cubPoints must have rank 3 (C,P,D).");
    INTREPID2_TEST_FOR_EXCEPTION( cubWeights.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolumeBoundary): cubWeights must have rank 2 (C,P).");
    INTREPID2_TEST_FOR_EXCEPTION( cellCoords.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolumeBoundary): cellCoords must have rank 3 of (C,P,D).");

    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.extent(0) != cellCoords.extent(0) ||
                                  cubPoints.extent(0) != cubWeights.extent(0), std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolume): cubPoints, cubWeights and cellCoords dimension(0) are not consistent, numCells");

    {
      const ordinal_type spaceDim = cellCoords.extent(2);
      const ordinal_type sideDim = spaceDim - 1;
      const size_type numPrimarySideNodes = primaryCellTopo_.getNodeCount(sideDim, sideIndex_);
      
      INTREPID2_TEST_FOR_EXCEPTION( cubPoints.extent(1) != numPrimarySideNodes || 
                                    cubWeights.extent(1) != numPrimarySideNodes, std::invalid_argument,
                                    ">>> ERROR (CubatureControlVolume): cubPoints and cubWeights dimension(1) are not consistent, numPrimarySideNodes");
    }
    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.extent(2) != cellCoords.extent(2) ||
                                  static_cast<ordinal_type>(cubPoints.extent(2)) != getDimension(), std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolume): cubPoints, cellCoords, this->getDimension() are not consistent, spaceDim.");
#endif
    
    typedef Kokkos::DynRankView<PT,SpT> tempPointViewType;
    typedef Kokkos::DynRankView<PT,Kokkos::LayoutStride,SpT> tempPointStrideViewType;

    // get array dimensions
    const ordinal_type numCells = cellCoords.extent(0);
    const ordinal_type numNodesPerCell = cellCoords.extent(1);
    const ordinal_type spaceDim = cellCoords.extent(2);
    const ordinal_type sideDim = spaceDim - 1;

    const ordinal_type numNodesPerSubcv = subcvCellTopo_.getNodeCount();
    tempPointViewType subcvCoords("CubatureControlVolumeBoundary::subcvCoords",
                                  numCells, numNodesPerCell, numNodesPerSubcv, spaceDim);
    CellTools<SpT>::getSubcvCoords(subcvCoords,
                                   cellCoords,
                                   primaryCellTopo_);

    //const auto numPrimarySides = primaryCellTopo_.getSubcellCount(sideDim);
    const ordinal_type numSubcvPoints = 1;

    tempPointViewType subcvJacobian("CubatureControlVolumeBoundary::subcvJacobian",
                                    numCells, numSubcvPoints, spaceDim, spaceDim);
    
    tempPointViewType subcvJacobianDet("CubatureControlVolumeBoundary::subcvJacobianDet",
                                       numCells, numSubcvPoints);
    
    tempPointViewType weights("CubatureControlVolumeBoundary::subcvWeights",
                              numCells, 1);
    Kokkos::deep_copy(weights, spaceDim == 2 ? 2.0 : 4.0);

    tempPointViewType scratch("CubatureControlVolumeBoundary::scratch",
                              numCells*numSubcvPoints*spaceDim*spaceDim);

    const ordinal_type numPrimarySideNodes = primaryCellTopo_.getNodeCount(sideDim, sideIndex_);
    for (ordinal_type node=0;node<numPrimarySideNodes;++node) {
      const auto sideRange = Kokkos::pair<ordinal_type,ordinal_type>(node, node+1);
      const auto sidePoint = Kokkos::subdynrankview(sidePoints_, sideRange, Kokkos::ALL());
      
      const auto idx = primaryCellTopo_.getNodeMap(sideDim, sideIndex_, node);
      auto subcvCoordsNode = Kokkos::subdynrankview(subcvCoords, Kokkos::ALL(), idx,  Kokkos::ALL(), Kokkos::ALL());
        
      CellTools<SpT>::setJacobian(subcvJacobian,        // C, P, D, D
                                  sidePoint,            //    P, D
                                  subcvCoordsNode,      // C, N, D
                                  subcvCellTopo_);
      
      CellTools<SpT>::setJacobianDet(subcvJacobianDet,  // C, P
                                     subcvJacobian);    // C, P, D, D
      
      auto cubPointsNode  = Kokkos::subdynrankview(cubPoints,  Kokkos::ALL(), node, Kokkos::ALL());
      
      typedef Kokkos::View<ordinal_type*,SpT> mapViewType;
      const auto sideNodeMap = Kokkos::subview(sideNodeMap_, node, Kokkos::ALL());
      
      typedef typename ExecSpace<typename pointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;
      
      const auto loopSize = numCells;
      Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
      
      // compute centers
      typedef Functor<pointViewType,tempPointStrideViewType,mapViewType> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(cubPointsNode,
                                                subcvCoordsNode,
                                                sideNodeMap) );
   
      // compute weights
      const auto sideOrd = boundarySidesHost_(sideIndex_, node);

      // cub weights node requires to have rank 2
      auto cubWeightsNode = Kokkos::subdynrankview(cubWeights, Kokkos::ALL(), sideRange);
      switch (spaceDim) {
      case 2: {
        std::cout << "subcv jacobian rank = " << subcvJacobian.rank() << std::endl;
        FunctionSpaceTools<SpT>::computeEdgeMeasure(cubWeightsNode, // rank 2
                                                    subcvJacobian,  // rank 4
                                                    weights,        // rank 2
                                                    sideOrd,
                                                    subcvCellTopo_,
                                                    scratch);
        break;
      }
      case 3: {
        FunctionSpaceTools<SpT>::computeFaceMeasure(cubWeightsNode,
                                                    subcvJacobian,
                                                    weights,
                                                    sideOrd,
                                                    subcvCellTopo_,
                                                    scratch);
        break;
      }
      }
    }
  }
} 

#endif

