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

namespace Intrepid {

template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureNodal<Scalar,ArrayPoint,ArrayWeight>::CubatureNodal(const shards::CellTopology & cellTopo) :
  cellTopo_(cellTopo.getBaseCellTopologyData()) {

  numPoints_ = cellTopo_.getVertexCount();
  cellDim_ = cellTopo_.getDimension();

  Scalar cellVolume(0);
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

template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureNodal<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint  & cubPoints,
                                                               ArrayWeight & cubWeights) const {
  // check size of cubPoints and cubWeights
  TEUCHOS_TEST_FOR_EXCEPTION( ( ( (unsigned)cubPoints.size() < numPoints_*cellDim_ ) || ( (unsigned)cubWeights.size() < numPoints_ ) ), std::out_of_range,
                              ">>> ERROR (CubatureNodal): Insufficient space allocated for cubature points or weights.");

  CellTools<Scalar>::getReferenceSubcellVertices(cubPoints, cellDim_, 0, cellTopo_);

  for (unsigned pointId = 0; pointId < numPoints_; pointId++) {
    cubWeights(pointId) = weightVal_;
  }

} // end getCubature

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureNodal<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint&  cubPoints,
                                                               ArrayWeight& cubWeights,
                                                               ArrayPoint&  cellCoords) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (CubatureNodal): Cubature defined in reference space calling method for physical space cubature.");
}



template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureNodal<Scalar,ArrayPoint,ArrayWeight>::getNumPoints() const {
  return numPoints_;
} // end getNumPoints



template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureNodal<Scalar,ArrayPoint,ArrayWeight>::getDimension() const {
  return cellDim_;
} // end dimension



template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureNodal<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(std::vector<int> & accuracy) const {
  accuracy.assign(1, 2);
} // end getAccuracy


} // end namespace Intrepid
