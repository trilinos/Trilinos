// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubatureControlVolumeDef.hpp
    \brief  Header file for the Intrepid2::CubatureControlVolume class.
    \author Created by K. Peterson, P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_CONTROLVOLUME_DEF_HPP__
#define __INTREPID2_CUBATURE_CONTROLVOLUME_DEF_HPP__

namespace Intrepid2 {

  template <typename DT, typename PT, typename WT>
  CubatureControlVolume<DT,PT,WT>::
  CubatureControlVolume(const shards::CellTopology cellTopology) {

    // define primary cell topology with given one
    primaryCellTopo_ = cellTopology;

    // subcell is defined either hex or quad according to dimension
    switch (primaryCellTopo_.getDimension()) {
    case 2:
      subcvCellTopo_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >());
      break;
    case 3:
      subcvCellTopo_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >());
      break;
    }

    // computation order is always one;
    degree_ = 1;

    // create subcell cubature points and weights and cache them
    const ordinal_type subcvDegree = 2;
    auto subcvCubature = DefaultCubatureFactory::create<DT,PT,WT>(subcvCellTopo_, subcvDegree);

    const ordinal_type numSubcvPoints = subcvCubature->getNumPoints();
    const ordinal_type subcvDim       = subcvCubature->getDimension();

    subcvCubaturePoints_  = Kokkos::DynRankView<PT,DT>("CubatureControlVolume::subcvCubaturePoints_",
                                                        numSubcvPoints, subcvDim);
    subcvCubatureWeights_ = Kokkos::DynRankView<WT,DT>("CubatureControlVolume::subcvCubatureWeights_",
                                                        numSubcvPoints);
    
    subcvCubature->getCubature(subcvCubaturePoints_,
                               subcvCubatureWeights_);
  }
  
  template <typename DT, typename PT, typename WT>
  void
  CubatureControlVolume<DT,PT,WT>::
  getCubature( PointViewType  cubPoints,
               weightViewType cubWeights,
               PointViewType  cellCoords ) const {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolume): cubPoints must have rank 3 (C,P,D).");
    INTREPID2_TEST_FOR_EXCEPTION( cubWeights.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolume): cubWeights must have rank 2 (C,P).");
    INTREPID2_TEST_FOR_EXCEPTION( cellCoords.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolume): cellCoords must have rank 3 of (C,P,D).");

    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.extent(0) != cellCoords.extent(0) ||
                                  cubPoints.extent(0) != cubWeights.extent(0), std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolume): cubPoints, cubWeights and cellCoords dimension(0) are not consistent, numCells");

    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.extent(1) != cellCoords.extent(1) ||
                                  cubPoints.extent(1) != cubWeights.extent(1), std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolume): cubPoints, cubWeights and cellCoords dimension(1) are not consistent, numNodesPerCell");

    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.extent(2) != cellCoords.extent(2) ||
                                  static_cast<ordinal_type>(cubPoints.extent(2)) != getDimension(), std::invalid_argument,
                                  ">>> ERROR (CubatureControlVolume): cubPoints, cellCoords, this->getDimension() are not consistent, spaceDim.");
#endif
    typedef Kokkos::DynRankView<PT,DT> tempPointViewType;

    // get array dimensions
    const ordinal_type numCells = cellCoords.extent(0);
    const ordinal_type numNodesPerCell = cellCoords.extent(1);
    const ordinal_type spaceDim = cellCoords.extent(2);

    const ordinal_type numNodesPerSubcv = subcvCellTopo_.getNodeCount();
    tempPointViewType subcvCoords("CubatureControlVolume::subcvCoords_",
                                  numCells, numNodesPerCell, numNodesPerSubcv, spaceDim);
    CellTools<DT>::getSubcvCoords(subcvCoords,
                                   cellCoords,
                                   primaryCellTopo_);

    const ordinal_type numSubcvPoints = subcvCubaturePoints_.extent(0);
    tempPointViewType subcvJacobian("CubatureControlVolume::subcvJacobian_",
                                    numCells, numNodesPerCell, numSubcvPoints, spaceDim, spaceDim);
    
    tempPointViewType subcvJacobianDet("CubatureControlVolume::subcvJacobDet_",
                                       numCells, numNodesPerCell, numSubcvPoints);
    
    // numNodesPerCell is maximum 8; this repeated run is necessary because of cell tools input consideration
    for (ordinal_type node=0;node<numNodesPerCell;++node) {
      auto subcvJacobianNode    = Kokkos::subdynrankview(subcvJacobian,    Kokkos::ALL(), node, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
      auto subcvCoordsNode      = Kokkos::subdynrankview(subcvCoords,      Kokkos::ALL(), node, Kokkos::ALL(), Kokkos::ALL());
      auto subcvJacobianDetNode = Kokkos::subdynrankview(subcvJacobianDet, Kokkos::ALL(), node, Kokkos::ALL());
      
      CellTools<DT>::setJacobian(subcvJacobianNode,    // C, P, D, D
                                  subcvCubaturePoints_, //    P, D
                                  subcvCoordsNode,      // C, N, D
                                  subcvCellTopo_);

      CellTools<DT>::setJacobianDet(subcvJacobianDetNode,  // C, P
                                     subcvJacobianNode);    // C, P, D, D
    }
    
    //typedef typename ExecSpace<typename PointViewType::execution_space,DT>::ExecSpaceType ExecSpaceType;
    
    const auto loopSize = numCells;
    Kokkos::RangePolicy<typename DT::execution_space,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    
    typedef Functor<PointViewType,weightViewType,tempPointViewType,tempPointViewType,tempPointViewType> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(cubPoints, 
                                              cubWeights, 
                                              subcvCoords, 
                                              subcvCubatureWeights_,
                                              subcvJacobianDet) );
  } 
    
} 

#endif

