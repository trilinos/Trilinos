// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubatureControlVolume.hpp
    \brief  Header file for the Intrepid2::CubatureControlVolume class.
    \author Created by K. Peterson, P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_CONTROLVOLUME_HPP__
#define __INTREPID2_CUBATURE_CONTROLVOLUME_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Cubature.hpp"

#include "Shards_CellTopology.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"

namespace Intrepid2{

  /** \class Intrepid2::CubatureControlVolume
      \brief Defines cubature (integration) rules over control volumes.

      Each primary cell contains one sub-control volume per node and
      there is one integration point per sub-control volume.
  */
  template<typename DeviceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubatureControlVolume
    : public Cubature<DeviceType,pointValueType,weightValueType> {
  public:

    template<typename cubPointViewType,
             typename cubWeightViewType,
             typename subcvCoordViewType,
             typename subcvWeightViewType,
             typename jacDetViewType>
    struct Functor {
      /**/  cubPointViewType    _cubPoints;
      /**/  cubWeightViewType   _cubWeights;
      const subcvCoordViewType  _subcvCoords;
      const subcvWeightViewType _subcvWeights;
      const jacDetViewType      _jacDets;

      KOKKOS_INLINE_FUNCTION
      Functor( cubPointViewType    cubPoints_,
               cubWeightViewType   cubWeights_,
               subcvCoordViewType  subcvCoords_,
               subcvWeightViewType subcvWeights_,
               jacDetViewType      jacDets_ )
        : _cubPoints(cubPoints_), _cubWeights(cubWeights_), 
          _subcvCoords(subcvCoords_), _subcvWeights(subcvWeights_), _jacDets(jacDets_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cell) const {        
        const ordinal_type numNodesPerCell  = _subcvCoords.extent(1);
        const ordinal_type numNodesPerSubcv = _subcvCoords.extent(2);
        const ordinal_type spaceDim         = _subcvCoords.extent(3);
        const ordinal_type numSubcvPoints   = _subcvWeights.extent(0);

        // compute subcv centers
        for (ordinal_type node=0;node<numNodesPerCell;++node) {
          typename cubPointViewType::value_type val[3] = {};
          for (ordinal_type subcv=0;subcv<numNodesPerSubcv;++subcv) {
            for (ordinal_type i=0;i<spaceDim;++i) 
              val[i] += _subcvCoords(cell, node, subcv, i);
          }
          for (ordinal_type i=0;i<spaceDim;++i) 
            _cubPoints(cell, node, i) = (val[i]/numNodesPerSubcv);
        }
        
        // compute weights (area or volume)
        for (ordinal_type node=0;node<numNodesPerCell;++node) {
          typename cubWeightViewType::value_type val = 0;
          for (ordinal_type pt=0;pt<numSubcvPoints;++pt)
            val += _subcvWeights(pt)*_jacDets(cell, node, pt);
          _cubWeights(cell, node) = val;
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
    Kokkos::DynRankView<pointValueType, DeviceType> subcvCubaturePoints_;
    Kokkos::DynRankView<weightValueType,DeviceType> subcvCubatureWeights_;

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
      return primaryCellTopo_.getNodeCount();
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
      return "CubatureControlVolume";
    }

    /** brief Constructor.
	\param cellTopology           [in]     - The topology of the primary cell.
    */
    CubatureControlVolume(const shards::CellTopology cellTopology);
    virtual ~CubatureControlVolume() {}

  }; // end class CubatureControlVolume

} // end namespace Intrepid2

#include "Intrepid2_CubatureControlVolumeDef.hpp"

#endif

