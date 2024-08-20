// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubatureControlVolumeBoundary.hpp
    \brief  Header file for the Intrepid2::CubatureControlVolumeBoundary class.
    \author Created by K. Peterson, P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_CONTROLVOLUME_BOUNDARY_HPP__
#define __INTREPID2_CUBATURE_CONTROLVOLUME_BOUNDARY_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Cubature.hpp"

#include "Shards_CellTopology.hpp"
#include "Intrepid2_CellTools.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"

namespace Intrepid2{
  
  /** \class Intrepid2::CubatureControlVolumeBoundary
      \brief Defines cubature (integration) rules over Neumann boundaries for control volume method.
      
      Integration on Neumann boundaries for the control volume method requires integration points
      defined on primary cell sides. These points are not equivalent to control volume points on lower
      dimensional topologies and therefore require a separate class to define them. 
  */
  template<typename DeviceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubatureControlVolumeBoundary
    : public Cubature<DeviceType,pointValueType,weightValueType> {
  public:

    template<typename cubPointViewType,
             typename subcvCoordViewType,
             typename mapViewType>
    struct Functor {
      /**/  cubPointViewType        _cubPoints;
      const subcvCoordViewType      _subcvCoords;
      const mapViewType             _sideMap;

      KOKKOS_INLINE_FUNCTION
      Functor( cubPointViewType        cubPoints_,
               subcvCoordViewType      subcvCoords_,
               mapViewType             sideMap_ )
        : _cubPoints(cubPoints_), 
          _subcvCoords(subcvCoords_), _sideMap(sideMap_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cell) const {        
        const ordinal_type numNodesPerSide = _sideMap(0);
        const ordinal_type spaceDim        = _cubPoints.extent(1);

        // compute side centers
        typename cubPointViewType::value_type val[3] = {};
        for (ordinal_type j=0;j<numNodesPerSide;++j) {
          for (ordinal_type i=0;i<spaceDim;++i) 
            val[i] += _subcvCoords(cell, _sideMap(j+1), i);
        }
        for (ordinal_type i=0;i<spaceDim;++i) 
          _cubPoints(cell, i) = (val[i]/numNodesPerSide);
        
      }
    };

  protected:
    
    /** \brief The topology of the primary cell side.
     */
    shards::CellTopology primaryCellTopo_;

    /** \brief The topology of the sub-control volume.
     */
    shards::CellTopology subcvCellTopo_;

    /** \brief The degree of the polynomials that are integrated exactly.
     */
    ordinal_type degree_;
    
    /** \brief Index of cell side
     */
    ordinal_type sideIndex_;

    // cubature points and weights associated with sub-control volume.
    Kokkos::View<ordinal_type**,Kokkos::HostSpace> boundarySidesHost_;
    Kokkos::View<ordinal_type**,Kokkos::LayoutRight,DeviceType> sideNodeMap_;
    Kokkos::DynRankView<pointValueType, DeviceType> sidePoints_;
    
  public:
    typedef typename Cubature<DeviceType,pointValueType,weightValueType>::PointViewType  PointViewType;
    typedef typename Cubature<DeviceType,pointValueType,weightValueType>::weightViewType weightViewType;

    using Cubature<DeviceType,pointValueType,weightValueType>::getCubature;

    /** \brief Returns cubature points and weights
        (return arrays must be pre-sized/pre-allocated).
	
	\param cubPoints             [out]        - Array containing the cubature points.
        \param cubWeights            [out]        - Array of corresponding cubature weights.
        \param cellCoords             [in]        - Array of cell coordinates
    */
    virtual
    void
    getCubature( PointViewType  cubPoints,
                 weightViewType cubWeights,
                 PointViewType  cellCoords) const override;
    
    /** \brief Returns the number of cubature points.
     */
    virtual
    ordinal_type
    getNumPoints() const override {
      // one control volume boundary cubature point per subcell node (for now)
      const ordinal_type sideDim = primaryCellTopo_.getDimension() - 1;
      return primaryCellTopo_.getNodeCount(sideDim, sideIndex_);
    }
    
    /** \brief Returns dimension of integration domain.
     */
    virtual
    ordinal_type 
    getDimension() const override {
      return primaryCellTopo_.getDimension();
    }
    
    /** \brief Returns cubature name.
     */
    virtual
    const char*
    getName() const override {
      return "CubatureControlVolumeBoundary";
    }

    /** brief Constructor.
	
	\param cellTopology           [in]     - The topology of the primary cell.
	\param cellSide               [in]     - The index of the boundary side of the primary cell 
    */
    CubatureControlVolumeBoundary(const shards::CellTopology cellTopology, 
                                  const ordinal_type         sideIndex);
    virtual ~CubatureControlVolumeBoundary() {}
    
  };

} 

#include "Intrepid2_CubatureControlVolumeBoundaryDef.hpp"

#endif

