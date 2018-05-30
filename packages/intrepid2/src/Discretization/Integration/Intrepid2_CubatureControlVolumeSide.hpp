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
  template<typename ExecSpaceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubatureControlVolumeSide
    : public Cubature<ExecSpaceType,pointValueType,weightValueType> {
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
    Kokkos::View<ordinal_type**,Kokkos::LayoutRight,ExecSpaceType> sideNodeMap_;
    Kokkos::DynRankView<pointValueType, ExecSpaceType> sidePoints_;

  public:
    typedef typename Cubature<ExecSpaceType,pointValueType,weightValueType>::pointViewType  pointViewType;
    typedef typename Cubature<ExecSpaceType,pointValueType,weightValueType>::weightViewType weightViewType;

    using Cubature<ExecSpaceType,pointValueType,weightValueType>::getCubature;

    /** \brief Returns cubature points and weights
        (return arrays must be pre-sized/pre-allocated).

	\param cubPoints             [out]        - Array containing the cubature points.
        \param cubWeights            [out]        - Array of corresponding cubature weights.
        \param cellCoords             [in]        - Array of cell coordinates
    */
    virtual
    void
    getCubature( pointViewType  cubPoints,
                 weightViewType cubWeights,
                 pointViewType  cellCoords) const;

    /** \brief Returns the number of cubature points.
     */
    virtual
    ordinal_type
    getNumPoints() const {
      return primaryCellTopo_.getEdgeCount();
    }

    /** \brief Returns dimension of integration domain.
     */
    virtual
    ordinal_type
    getDimension() const {
      return primaryCellTopo_.getDimension();
    }

    /** \brief Returns cubature name.
     */
    virtual
    const char*
    getName() const {
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

