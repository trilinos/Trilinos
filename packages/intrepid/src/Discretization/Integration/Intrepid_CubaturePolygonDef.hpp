// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubaturePolygonDef.hpp
    \brief  Definition file for the Intrepid::CubaturePolygon class.
    \author Created by P. Bochev and J. Lai.
*/


#include "Intrepid_CubatureDirectTriDefault.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include <vector>
#include <iostream>

namespace Intrepid{

  template<class Scalar, class ArrayPoint, class ArrayWeight>
  CubaturePolygon<Scalar,ArrayPoint,ArrayWeight>::CubaturePolygon(const shards::CellTopology& cellTopology,
								  const ArrayPoint& cellVertices,
								  int degree)
    : degree_(degree), cubDimension_(2), cellTopology_(cellTopology), cellVertices_(cellVertices){
    
    TEUCHOS_TEST_FOR_EXCEPTION( (degree < 0) || degree > INTREPID_CUBATURE_TRI_DEFAULT_MAX, std::out_of_range,
			">>> ERROR (CubaturePolygon): No direct cubature rule implemented for the desired polynomial degree.");
    // compute area and centroid of polygon
    Scalar area;
    std::vector<Scalar> centroid(2,0);
    int numNodes = cellTopology_.getNodeCount();
    
    for (int i=0;i<numNodes;i++){
      int first = cellTopology_.getNodeMap(1,i,0);
      int second = cellTopology_.getNodeMap(1,i,1);
      area += cellVertices_(first,0)*cellVertices_(second,1) - cellVertices_(second,0)*cellVertices_(first,1);
      centroid[0] += (cellVertices_(first,0) + cellVertices_(second,0))*(cellVertices_(first,0)*cellVertices_(second,1) - cellVertices_(second,0)*cellVertices_(first,1));
      centroid[1] += (cellVertices_(first,1) + cellVertices_(second,1))*(cellVertices_(first,0)*cellVertices_(second,1) - cellVertices_(second,0)*cellVertices_(first,1));
    }
    area /= 2;
    centroid[0] /= (6*area);
    centroid[1] /= (6*area);
        
    // get cubature for reference triangle
    CubatureDirectTriDefault<Scalar,ArrayPoint,ArrayWeight> cubatureTri(degree_);
    int numCubPointsPerTri = cubatureTri.getNumPoints();
    int cubDim = cubatureTri.getDimension();
    cubDimension_ = cubDim;
    FieldContainer<Scalar> cubatureTriPoints(numCubPointsPerTri,cubDim);
    
    FieldContainer<Scalar> cubatureTriWeights(numCubPointsPerTri);
    cubatureTri.getCubature(cubatureTriPoints,cubatureTriWeights);
    
    // copy into (C,P,D) sized field container where C is the number of triangles in polygon
    int numCells = cellTopology_.getEdgeCount();
    FieldContainer<Scalar> cubatureCellPoints(numCells,numCubPointsPerTri,cubDim);
    for (int k=0;k<numCells;k++){
      for (int i=0;i<numCubPointsPerTri;i++){
	for (int j=0;j<cubDim;j++){
	  cubatureCellPoints(k,i,j) = cubatureTriPoints(i,j);
	}
      }
    }
    
    
    // now map cubature to each triangle cell
    shards::CellTopology triangleTopology(shards::getCellTopologyData<shards::Triangle<3> >());
    int totalCubPoints = numCubPointsPerTri*cellTopology_.getEdgeCount();
    numPoints_ = totalCubPoints;
    cubaturePoints_.resize(totalCubPoints,cubDim);
    cubatureWeights_.resize(totalCubPoints);
    
    FieldContainer<Scalar> physicalPoints(numCells,numCubPointsPerTri,cubDim);
    FieldContainer<Scalar> trianglePoints(numCells,3,cubDim);
    int currPoint = 0;
    for (int i=0;i<numCells;i++){
      for (int j=0;j<cubDim;j++){
	trianglePoints(i,0,j) = cellVertices_(cellTopology_.getNodeMap(1,i,0),j);
	trianglePoints(i,1,j) = cellVertices_(cellTopology_.getNodeMap(1,i,1),j);
	trianglePoints(i,2,j) = centroid[j];
      }
    }
    
    CellTools<Scalar>::mapToPhysicalFrame(physicalPoints,cubatureTriPoints,trianglePoints,triangleTopology);
    
    // compute area of each triangle cell -- need when computing new weights
    FieldContainer<Scalar> jacobians(numCells,numCubPointsPerTri,cubDim,cubDim);
    FieldContainer<Scalar> detJacobians(numCells, numCubPointsPerTri);

    CellTools<Scalar>::setJacobian(jacobians,physicalPoints,trianglePoints,triangleTopology);
    CellTools<Scalar>::setJacobianDet(detJacobians,jacobians);
    
    for (int i=0;i<numCells;i++){
      for (int j=0;j<numCubPointsPerTri;j++){
	for (int k=0;k<cubDim;k++){
	  cubaturePoints_(currPoint,k) = physicalPoints(i,j,k);
	}
	cubatureWeights_(currPoint++) = cubatureTriWeights(j)*detJacobians(i,j);
      }
    }
  }  // end Constructor


template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubaturePolygon<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
								 ArrayWeight& cubWeights)const{
  int numCubPoints = numPoints_;
  int cellDim = cubDimension_;
    
  TEUCHOS_TEST_FOR_EXCEPTION ( ( cubPoints.size() < numCubPoints*cellDim || cubWeights.size() < numCubPoints ),
		       std::out_of_range,
		       ">>> ERROR (CubaturePolygon): Insufficient space allocated for cubature points or weights.");

  for (int pointId = 0; pointId < numCubPoints; pointId++){
    for (int dim = 0; dim < cellDim; dim++){
      cubPoints(pointId,dim) = cubaturePoints_(pointId,dim);
    }
    cubWeights(pointId) = cubatureWeights_(pointId);
  }
} // end getCubature

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubaturePolygon<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
                                                                 ArrayWeight& cubWeights,
                                                                 ArrayPoint& cellCoords) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (CubaturePolygon): Cubature defined in reference space calling method for physical space cubature.");
}

  

template<class Scalar, class ArrayPoint, class ArrayWeight>
int CubaturePolygon<Scalar,ArrayPoint,ArrayWeight>::getNumPoints()const{
  return numPoints_;
} // end getNumPoints

template<class Scalar, class ArrayPoint, class ArrayWeight>
int CubaturePolygon<Scalar,ArrayPoint,ArrayWeight>::getDimension()const{
  return cubDimension_;
} // end getDimension

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubaturePolygon<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(std::vector<int>& accuracy)const{
  accuracy.assign(1,degree_);
} // end getAccuracy
  
} // namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

