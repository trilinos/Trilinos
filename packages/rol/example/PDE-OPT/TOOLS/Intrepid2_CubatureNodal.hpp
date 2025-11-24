// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
//
//                           Intrepid Package
// Copyright 2007 NTESS and the Intrepid contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid_CubatureNodal.hpp
    \brief  Header file for the Intrepid::CubatureNodal class.
    \author Created by D. Ridzal.
*/

#ifndef INTREPID2_CUBATURE_NODAL_HPP
#define INTREPID2_CUBATURE_NODAL_HPP

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_Cubature.hpp"
#include "Teuchos_Assert.hpp"


namespace Intrepid2 {

/** \class Intrepid::CubatureNodal
    \brief Defines nodal cubature on select cell topologies, often used for
           lower-order integration or underintegration.  The use of this
           cubature with P1 elements yields a diagonal mass matrix.
*/
template<typename DeviceType=void, typename PointValueType=double, typename WeightValueType=double>
class CubatureNodal : public Cubature<DeviceType,PointValueType,WeightValueType> {
  public:
    using ExecSpaceType = typename DeviceType::execution_space;
    using PointViewType             = Kokkos::DynRankView<PointValueType,Kokkos::LayoutStride,DeviceType>;
    using WeightViewType            = Kokkos::DynRankView<WeightValueType,Kokkos::LayoutStride,DeviceType>;

    using PointViewTypeAllocatable  = Kokkos::DynRankView<PointValueType,DeviceType>;  // uses default layout; allows us to allocate (in contrast to LayoutStride)
    using WeightViewTypeAllocatable = Kokkos::DynRankView<WeightValueType,DeviceType>; // uses default layout; allows us to allocate (in contrast to LayoutStride)
    using TensorPointDataType       = TensorPoints<PointValueType,DeviceType>;
    using TensorWeightDataType      = TensorData<WeightValueType,DeviceType>;
  private:
  /** \brief  Base topology of the cells for which the cubature is defined. See  the Shards package
              http://trilinos.sandia.gov/packages/shards for definition of base cell topology.
  */
  shards::CellTopology cellTopo_;

  unsigned numPoints_;
  unsigned cellDim_;

  WeightValueType weightVal_;

  public:

  ~CubatureNodal() {}

  /** \brief Constructor.
       \param cellTopo             [in]        - Cell topology.
  */
  CubatureNodal(const shards::CellTopology & cellTopo);

  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

      \param cubPoints       [out]     - Vector containing the cubature points.
      \param cubWeights      [out]     - Vector of corresponding cubature weights.
  */
  virtual void getCubature(PointViewType cubPoints,
                           WeightViewType cubWeights) const override;

  /** \brief Returns cubature points and weights.
              Method for physical space cubature, throws an exception.

       \param cubPoints             [out]        - Array containing the cubature points.
       \param cubWeights            [out]        - Array of corresponding cubature weights.
       \param cellCoords             [in]        - Array of cell coordinates
  */
  virtual void getCubature(PointViewType cubPoints,
                           WeightViewType cubWeights,
                           PointViewType cellCoords) const override;

  /** \brief Returns the number of cubature points.
  */
  virtual ordinal_type getNumPoints() const override;

  /** \brief Returns dimension of integration domain.
  */
  virtual ordinal_type getDimension() const override;

  /** \brief Returns max. degree of polynomials that are integrated exactly.
             The return vector has the size of the degree_ vector.
  */
  virtual ordinal_type getAccuracy() const override;

}; // end class CubatureNodal


} // end namespace Intrepid


// include templated definitions
#include "Intrepid2_CubatureNodalDef.hpp"

#endif
