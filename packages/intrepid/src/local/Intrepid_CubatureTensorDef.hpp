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

/** \file   Intrepid_CubatureTensorDef.hpp
    \brief  Definition file for the Intrepid::CubatureTensor class.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

template <class Scalar>
CubatureTensor<Scalar>::CubatureTensor(const ECell                      cellType,
                                       const int                        degree) :
  cellType_(cellType), degree_(degree) {
  switch (cellType_) {

    case CELL_QUAD:
      TEST_FOR_EXCEPTION((degree_ < 0) || (degree_ > INTREPID_MAX_CUBATURE_DEGREE_EDGE),
                         std::out_of_range,
                         ">>> ERROR (CubatureTensor): No tensor-product cubature rule implemented for the desired polynomial degree.");
      break;

    case CELL_HEX:
      TEST_FOR_EXCEPTION((degree_ < 0) || (degree_ > INTREPID_MAX_CUBATURE_DEGREE_EDGE),
                         std::out_of_range,
                         ">>> ERROR (CubatureTensor): No tensor-product cubature rule implemented for the desired polynomial degree.");
      break;

    case CELL_TRIPRISM:
      TEST_FOR_EXCEPTION((degree_ < 0) || (degree_ > std::min(INTREPID_MAX_CUBATURE_DEGREE_EDGE,INTREPID_MAX_CUBATURE_DEGREE_TRI)),
                         std::out_of_range,
                         ">>> ERROR (CubatureTensor): No tensor-product cubature rule implemented for the desired polynomial degree.");
      break;

    default:
       TEST_FOR_EXCEPTION((cellType_ != CELL_QUAD) && (cellType_ != CELL_HEX) && (cellType_ != CELL_TRIPRISM),
                          std::invalid_argument,
                          ">>> ERROR (CubatureTensor): Invalid cell type.");
  } // end switch
}



template <class Scalar>
void CubatureTensor<Scalar>::getCubature(int &                            numCubPoints,
                                         Teuchos::Array< Point<Scalar> >& cubPoints,
                                         Teuchos::Array<Scalar>&          cubWeights) const {

  numCubPoints = getNumPoints();
  
  int cellDim = MultiCell<Scalar>::getCellDim(cellType_);

  Point<Scalar> tempPoint(cellDim);
  cubPoints.assign(numCubPoints,tempPoint);
  cubWeights.assign(numCubPoints,(Scalar)0);

  getCubature(cubPoints, cubWeights);
} // end getCubature



template <class Scalar>
void CubatureTensor<Scalar>::getCubature(Teuchos::Array< Point<Scalar> >& cubPoints,
                                         Teuchos::Array<Scalar>&          cubWeights) const {

  EFrame cubPointFrame = FRAME_REFERENCE;

  int numCubPoints = getNumPoints();

  TEST_FOR_EXCEPTION( ( ( (int)cubPoints.size() < numCubPoints ) || ( (int)cubWeights.size() < numCubPoints ) ),
                        std::out_of_range,
                        ">>> ERROR (CubatureTensor): Insufficient space allocated for cubature points or weights.");

  Scalar x[3];
  int d[3];
 
  switch (cellType_) {

    case CELL_QUAD: {
        //
        // Find the digits d[0]-d[dim-1] of the point_id as a base numEdgePoints
        // number. The value of digit d[i] gives the 1D cubature point to use for
        // the i-th spatial coordinate of the tensor product cubature point.
        // Cubature weights are computed as products of weights of lower-dimensional
        // rules.
        //
        int cubatureIndexEdge = CUBATURE_GAUSS_0 + degree_;
        int numEdgePoints     = (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].numPoints_;
        for(int point_id = 0; point_id < numCubPoints; point_id++){
          d[0] = point_id / numEdgePoints;
          d[1] = point_id % numEdgePoints;
          cubWeights[point_id] = 1.0;
          for(int dim = 0; dim < 2; dim++) {
            x[dim]                = (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].points_[d[dim]][0];
            cubWeights[point_id] *= (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].weights_[d[dim]];
          }
          cubPoints[point_id].setCoordinates(x,2);
          cubPoints[point_id].setFrameKind(cubPointFrame);
        }
      }
      break;

    case CELL_HEX: {
        //
        // Find the digits d[0]-d[dim-1] of the point_id as a base numEdgePoints
        // number. The value of digit d[i] gives the 1D cubature point to use for
        // the i-th spatial coordinate of the tensor product cubature point.
        // Cubature weights are computed as products of weights of lower-dimensional
        // rules.
        //
        int cubatureIndexEdge = CUBATURE_GAUSS_0 + degree_;
        int numEdgePoints     = (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].numPoints_;
        int numEdgePointsSq   = numEdgePoints*numEdgePoints;
        int point_id_reduced  = 0;
        for(int point_id = 0; point_id < numCubPoints; point_id++){
          point_id_reduced = point_id % numEdgePointsSq;
          d[0] = point_id / numEdgePointsSq;
          d[1] = point_id_reduced / numEdgePoints;
          d[2] = point_id_reduced % numEdgePoints;
          cubWeights[point_id] = 1.0;
          for(int dim = 0; dim < 3; dim++) {
            x[dim]                = (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].points_[d[dim]][0];
            cubWeights[point_id] *= (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].weights_[d[dim]];
          }
          cubPoints[point_id].setCoordinates(x,3);
          cubPoints[point_id].setFrameKind(cubPointFrame);
        }
      }
      break;

    case CELL_TRIPRISM: {
        //
        // Find the digits d[0]-d[dim-1] of the point_id as a base numTriPoints,numEdgePoints
        // number. The value of digit d[1] gives the 1D Gauss point to use for
        // the 1st spatial coordinate of the tensor product cubature point.
        // The value of digit d[0] gives the 2D Triangle point to use for the
        // 2nd and 3rd spatial coordinates.
        // Cubature weights are computed as products of weights of lower-dimensional
        // rules.
        //
        int cubatureIndexEdge = CUBATURE_GAUSS_0 + degree_;
        int cubatureIndexTri  = CUBATURE_TRI_0   + degree_;
        //int numEdgePoints     = (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].numPoints_;
        int numTriPoints      = (CubatureDirect<Scalar>::exposeData())[cubatureIndexTri].numPoints_;
        for(int point_id = 0; point_id < numCubPoints; point_id++){
          d[0] = point_id / numTriPoints;
          d[1] = point_id % numTriPoints;
          x[0]                  = (CubatureDirect<Scalar>::exposeData())[cubatureIndexTri].points_[d[1]][0];
          x[1]                  = (CubatureDirect<Scalar>::exposeData())[cubatureIndexTri].points_[d[1]][1];
          x[2]                  = (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].points_[d[0]][0];
          cubWeights[point_id]  = (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].weights_[d[0]] *
                                  (CubatureDirect<Scalar>::exposeData())[cubatureIndexTri].weights_[d[1]];
          cubPoints[point_id].setCoordinates(x,3);
          cubPoints[point_id].setFrameKind(cubPointFrame);
        }
      }
      break;

    default:
       TEST_FOR_EXCEPTION((cellType_ != CELL_QUAD) && (cellType_ != CELL_HEX) && (cellType_ != CELL_TRIPRISM),
                          std::invalid_argument,
                          ">>> ERROR (CubatureTensor): Invalid cell type.");
  } // end switch
} // end getCubature



template<class Scalar>
int CubatureTensor<Scalar>::getNumPoints() const {

  int numCubPoints = -1;

  switch (cellType_) {

    case CELL_QUAD: {
        int cubatureIndexEdge = CUBATURE_GAUSS_0 + degree_;
        int numEdgePoints     = (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].numPoints_;
        numCubPoints          = numEdgePoints*numEdgePoints;
      }
      break;

    case CELL_HEX: {
        int cubatureIndexEdge = CUBATURE_GAUSS_0 + degree_;
        int numEdgePoints     = (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].numPoints_;
        numCubPoints          = numEdgePoints*numEdgePoints*numEdgePoints;
      }
      break;

    case CELL_TRIPRISM: {
        int cubatureIndexEdge = CUBATURE_GAUSS_0 + degree_;
        int cubatureIndexTri  = CUBATURE_TRI_0   + degree_;
        int numEdgePoints     = (CubatureDirect<Scalar>::exposeData())[cubatureIndexEdge].numPoints_;
        int numTriPoints      = (CubatureDirect<Scalar>::exposeData())[cubatureIndexTri].numPoints_;
        numCubPoints          = numEdgePoints*numTriPoints;
      }
      break;

    default:
       TEST_FOR_EXCEPTION((cellType_ != CELL_QUAD) && (cellType_ != CELL_HEX) && (cellType_ != CELL_TRIPRISM),
                          std::invalid_argument,
                          ">>> ERROR (CubatureTensor): Invalid cell type.");
  } // end switch

  return numCubPoints;

} // end getNumPoints



template <class Scalar>
ECell CubatureTensor<Scalar>::getCellType() const {
  return cellType_;
}



template <class Scalar>
int CubatureTensor<Scalar>::getAccuracy() const {
  return degree_;
}

} // end namespace Intrepid
