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

/** \file   Intrepid_CubatureNodalDef.hpp
    \brief  Definition file for the Intrepid::CubatureNodal class.
    \author Created by D. Ridzal.
*/

namespace Intrepid2 {

template <class DeviceType, class PointValueType, class WeightValueType>
CubatureNodal<DeviceType,PointValueType,WeightValueType>::CubatureNodal(const shards::CellTopology & cellTopo) :
  cellTopo_(cellTopo.getBaseCellTopologyData()) {

  numPoints_ = cellTopo_.getVertexCount();
  cellDim_ = cellTopo_.getDimension();

  WeightValueType cellVolume(0);
  switch( cellTopo_.getKey() ){
    case shards::Line<2>::key:
      cellVolume = 2.0;  break;
    case shards::Triangle<3>::key:
      cellVolume = 0.5;  break;
    case shards::Quadrilateral<4>::key:
      cellVolume = 4.0;  break;
    case shards::Tetrahedron<4>::key:
      cellVolume = 1.0/6.0;  break;
    case shards::Hexahedron<8>::key:
      cellVolume = 8.0;  break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( (true), std::invalid_argument,
                                 ">>> ERROR (Intrepid::CubatureNodal): Cell topology not supported.");
  } // switch key

  weightVal_ = cellVolume/numPoints_;
}

template <class DeviceType, class PointValueType, class WeightValueType>
void CubatureNodal<DeviceType,PointValueType,WeightValueType>::getCubature(PointViewType  cubPoints,
                                                                           WeightViewType cubWeights) const {
  // check size of cubPoints and cubWeights
  TEUCHOS_TEST_FOR_EXCEPTION( ( ( (unsigned)cubPoints.size() < numPoints_*cellDim_ ) || ( (unsigned)cubWeights.size() < numPoints_ ) ), std::out_of_range,
                              ">>> ERROR (CubatureNodal): Insufficient space allocated for cubature points or weights.");

  CellTools<DeviceType>::getReferenceSubcellVertices(cubPoints, cellDim_, 0, cellTopo_);

  for (unsigned pointId = 0; pointId < numPoints_; pointId++) {
    cubWeights(pointId) = weightVal_;
  }

} // end getCubature

template <class DeviceType, class PointValueType, class WeightValueType>
void CubatureNodal<DeviceType,PointValueType,WeightValueType>::getCubature(PointViewType  cubPoints,
                                                                           WeightViewType cubWeights,
                                                                           PointViewType  cellCoords) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (CubatureNodal): Cubature defined in reference space calling method for physical space cubature.");
}



template <class DeviceType, class PointValueType, class WeightValueType>
ordinal_type CubatureNodal<DeviceType,PointValueType,WeightValueType>::getNumPoints() const {
  return numPoints_;
} // end getNumPoints



template <class DeviceType, class PointValueType, class WeightValueType>
ordinal_type CubatureNodal<DeviceType,PointValueType,WeightValueType>::getDimension() const {
  return cellDim_;
} // end dimension



template <class DeviceType, class PointValueType, class WeightValueType>
ordinal_type CubatureNodal<DeviceType,PointValueType,WeightValueType>::getAccuracy() const {
  return 2;
} // end getAccuracy


} // end namespace Intrepid
