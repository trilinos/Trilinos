// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubatureControlVolumeSide.hpp
    \brief  Header file for the Intrepid2::CubatureControlVolumeSide class.
    \author Created by K. Peterson, P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_CONTROLVOLUME_SIDE_HPP__
#define __INTREPID2_CUBATURE_CONTROLVOLUME_SIDE_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Cubature.hpp"

#include "Shards_CellTopology.hpp"
#include "Intrepid2_CellTools.hpp"

namespace Intrepid2{

  /** \class Intrepid2::CubatureControlVolumeSide
      \brief Defines cubature (integration) rules over control volumes.
  */
  template<typename DeviceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubatureControlVolumeSide
    : public Cubature<DeviceType,pointValueType,weightValueType> {
  public:

    template<typename cubPointViewType,
             typename cubWeightViewType,
             typename subcvCoordViewType,
             typename subcvSideNormalViewType,
             typename mapViewType>
    struct Functor {
            cubPointViewType        _cubPoints;
            cubWeightViewType       _cubWeights;
      const subcvCoordViewType      _subcvCoords;
      const subcvSideNormalViewType _subcvSideNormals;
      const mapViewType             _sideMap;

      KOKKOS_INLINE_FUNCTION
      Functor( cubPointViewType        cubPoints_,
               cubWeightViewType       cubWeights_,
               subcvCoordViewType      subcvCoords_,
               subcvSideNormalViewType subcvSideNormals_,
               mapViewType             sideMap_ )
        : _cubPoints(cubPoints_), _cubWeights(cubWeights_), 
          _subcvCoords(subcvCoords_), _subcvSideNormals(subcvSideNormals_), _sideMap(sideMap_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cell) const {        
        const ordinal_type numNodesPerCell  = _cubPoints.extent(1);
        const ordinal_type spaceDim         = _cubPoints.extent(2);

        const ordinal_type numNodesPerSide  = _sideMap(0);
        const ordinal_type numSubcvPoints   = _subcvSideNormals.extent(2);

        const ordinal_type sideDim = spaceDim - 1;

        // compute side centers
        for (ordinal_type node=0;node<numNodesPerCell;++node) {
          typename cubPointViewType::value_type val[3] = {};
          for (ordinal_type j=0;j<numNodesPerSide;++j) {
            for (ordinal_type i=0;i<spaceDim;++i) 
              val[i] += _subcvCoords(cell, node, _sideMap(j+1), i);
          }
          for (ordinal_type i=0;i<spaceDim;++i) 
            _cubPoints(cell, node, i) = (val[i]/numNodesPerSide);
        }
        
        // compute weights (area or volume)
        for (ordinal_type node=0;node<numNodesPerCell;++node) {
          for (ordinal_type i=0;i<spaceDim;++i) {
            typename cubWeightViewType::value_type val = 0;
            for (ordinal_type pt=0;pt<numSubcvPoints;++pt)
              val += _subcvSideNormals(cell, node, pt, i)*pow(2,sideDim);
            _cubWeights(cell, node, i) = val;
          }
        }
      }
    };

  protected:

    /** \brief The topology of the primary cell.
     */
    shards::CellTopology primaryCellTopo_;

    /** \brief The topology of the sub-control volume.
     */
    shards::CellTopology subcvCellTopo_;

    /** \brief The degree of the polynomials that are integrated exactly.
     */
    ordinal_type degree_;

    // cubature points and weights associated with sub-control volume.
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
      return primaryCellTopo_.getEdgeCount();
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
      return "CubatureControlVolumeSide";
    }

    /** brief Constructor.
	\param cellTopology           [in]     - The topology of the primary cell.
    */
    CubatureControlVolumeSide(const shards::CellTopology cellTopology);
    virtual ~CubatureControlVolumeSide() {}

  }; // end class CubatureControlVolume

} // end namespace Intrepid2

#include "Intrepid2_CubatureControlVolumeSideDef.hpp"

#endif

