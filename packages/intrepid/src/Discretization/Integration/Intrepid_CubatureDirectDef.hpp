// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubatureDirectDef.hpp
    \brief  Definition file for the Intrepid::CubatureDirect class.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureDirect<Scalar,ArrayPoint,ArrayWeight>::getCubatureData(ArrayPoint  &                cubPoints,
                                                                    ArrayWeight &                cubWeights,
                                                                    const CubatureTemplate *     cubData) const {

  int numCubPoints = getNumPoints();
  int cellDim      = getDimension();
  // check size of cubPoints and cubWeights
  TEST_FOR_EXCEPTION( ( ( (int)cubPoints.size() < numCubPoints*cellDim ) || ( (int)cubWeights.size() < numCubPoints ) ),
                      std::out_of_range,
                      ">>> ERROR (CubatureDirect): Insufficient space allocated for cubature points or weights.");

  for (int pointId = 0; pointId < numCubPoints; pointId++) {
    for (int dim = 0; dim < cellDim; dim++) {
      cubPoints(pointId,dim) = cubData->points_[pointId][dim];
    }
    cubWeights(pointId) = cubData->weights_[pointId];
  }
} // end getCubatureData



template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureDirect<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint  & cubPoints,
                                                                ArrayWeight & cubWeights) const {
  getCubatureData( cubPoints, cubWeights, &(exposeCubatureData()[degree_]) );
} // end getCubature



template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureDirect<Scalar,ArrayPoint,ArrayWeight>::getNumPoints() const {
  return exposeCubatureData()[degree_].numPoints_;
} // end getNumPoints



template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureDirect<Scalar,ArrayPoint,ArrayWeight>::getDimension() const {
  return dimension_;
} // end dimension



template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureDirect<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(std::vector<int> & accuracy) const {
  accuracy.assign(1, degree_);
} // end getAccuracy

    
} // end namespace Intrepid
