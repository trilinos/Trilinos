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

/** \file   Intrepid_CubatureTensorDef.hpp
    \brief  Definition file for the Intrepid::CubatureTensor class.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::CubatureTensor(std::vector< Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > > cubatures) {
  unsigned numCubs = cubatures.size();
  TEUCHOS_TEST_FOR_EXCEPTION( (numCubs < 1),
                      std::out_of_range,
                      ">>> ERROR (CubatureTensor): Input cubature array must be of size 1 or larger.");

  cubatures_ = cubatures;

  unsigned numDegrees = 0;
  for (unsigned i=0; i<numCubs; i++) {
    std::vector<int> tmp;
    cubatures[i]->getAccuracy(tmp);
    numDegrees += tmp.size();
  }

  degree_.assign(numDegrees, 0);
  int count  = 0;
  dimension_ = 0;
  for (unsigned i=0; i<numCubs; i++) {
    std::vector<int> tmp;
    cubatures[i]->getAccuracy(tmp);
    for (unsigned j=0; j<tmp.size(); j++) {
      degree_[count] = tmp[j];
      count++;
    }
    dimension_ += cubatures[i]->getDimension();
  }
}



template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::CubatureTensor(Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature1,
                                                              Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature2) {
  cubatures_.resize(2);
  cubatures_[0] = cubature1;
  cubatures_[1] = cubature2;

  degree_.assign(2, 0);
  for (unsigned i=0; i<2; i++){
      std::vector<int> d;
      cubatures_[i]->getAccuracy(d); degree_[i] = d[0];
  }

  dimension_ = cubatures_[0]->getDimension() + cubatures_[1]->getDimension();
}



template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::CubatureTensor(Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature1,
                                                              Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature2,
                                                              Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature3) {
  cubatures_.resize(3);
  cubatures_[0] = cubature1;
  cubatures_[1] = cubature2;
  cubatures_[2] = cubature3;

  degree_.assign(3, 0);
  for (unsigned i=0; i<3; i++){
      std::vector<int> d;
      cubatures_[i]->getAccuracy(d); degree_[i] = d[0];
  }

  dimension_ = cubatures_[0]->getDimension() + cubatures_[1]->getDimension() + cubatures_[2]->getDimension();
}



template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::CubatureTensor(Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature, int n) {
  cubatures_.resize(n);
  for (int i=0; i<n; i++) {
    cubatures_[i] = cubature;
  }

  std::vector<int> d;
  cubatures_[0]->getAccuracy(d);
  degree_.assign(n,d[0]);

  dimension_ = cubatures_[0]->getDimension()*n;
}



template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint  & cubPoints,
                                                                ArrayWeight & cubWeights) const {
  int numCubPoints = getNumPoints();
  int cubDim       = getDimension();
  // check size of cubPoints and cubWeights
  TEUCHOS_TEST_FOR_EXCEPTION( ( ( (int)cubPoints.size() < numCubPoints*cubDim ) || ( (int)cubWeights.size() < numCubPoints ) ),
                      std::out_of_range,
                      ">>> ERROR (CubatureTensor): Insufficient space allocated for cubature points or weights.");

  unsigned numCubs   = cubatures_.size();
  std::vector<unsigned> numLocPoints(numCubs);
  std::vector<unsigned> locDim(numCubs);
  std::vector< FieldContainer<Scalar> > points(numCubs);
  std::vector< FieldContainer<Scalar> > weights(numCubs);

  // extract required points and weights
  for (unsigned i=0; i<numCubs; i++) {

    numLocPoints[i] = cubatures_[i]->getNumPoints();
    locDim[i]       = cubatures_[i]->getDimension();
    points[i].resize(numLocPoints[i], locDim[i]);
    weights[i].resize(numLocPoints[i]);

    // cubPoints and cubWeights are used here only for temporary data retrieval
    cubatures_[i]->getCubature(cubPoints, cubWeights);
    for (unsigned pt=0; pt<numLocPoints[i]; pt++) {
      for (unsigned d=0; d<locDim[i]; d++) {
        points[i](pt,d) = cubPoints(pt,d);
        weights[i](pt)  = cubWeights(pt);
      }
    }

  }

  // reset all weights to 1.0
  for (int i=0; i<numCubPoints; i++) {
      cubWeights(i) = (Scalar)1.0;
  }

  // fill tensor-product cubature
  int globDimCounter = 0;
  int shift          = 1;
  for (unsigned i=0; i<numCubs; i++) {

    for (int j=0; j<numCubPoints; j++) {
      /* int itmp = ((j*shift) % numCubPoints) + (j / (numCubPoints/shift)); // equivalent, but numerically unstable */
      int itmp = (j % (numCubPoints/shift))*shift + (j / (numCubPoints/shift));
      for (unsigned k=0; k<locDim[i]; k++) {
        cubPoints(itmp , globDimCounter+k) = points[i](j % numLocPoints[i], k);
      }
      cubWeights( itmp ) *= weights[i](j % numLocPoints[i]);
    }
    
    shift *= numLocPoints[i];
    globDimCounter += locDim[i];
  }

} // end getCubature

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
                                                                ArrayWeight& cubWeights,
                                                                ArrayPoint& cellCoords) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (CubatureTensor): Cubature defined in reference space calling method for physical space cubature.");
}


template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::getNumPoints() const {
  unsigned numCubs = cubatures_.size();
  int numCubPoints = 1;
  for (unsigned i=0; i<numCubs; i++) {
    numCubPoints *= cubatures_[i]->getNumPoints();
  }
  return numCubPoints;
} // end getNumPoints


template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::getDimension() const {
  return dimension_;
} // end dimension



template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(std::vector<int> & degree) const {
  degree = degree_;
} // end getAccuracy

} // end namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

