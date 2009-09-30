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

template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::CubatureTensor(std::vector< Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > > cubatures) {
  unsigned numCubs = cubatures.size();
  TEST_FOR_EXCEPTION( (numCubs < 1),
                      std::out_of_range,
                      ">>> ERROR (CubatureTensor): Input cubature array must be of size 1 or larger.");

  cubatures_ = cubatures;

  std::vector<int> tmp;
  unsigned numDegrees = 0;
  for (unsigned i=0; i<numCubs; i++) {
    cubatures[i]->getAccuracy(tmp);
    numDegrees += tmp.size();
  }

  degree_.assign(numDegrees, 0);
  int count  = 0;
  dimension_ = 0;
  for (unsigned i=0; i<numCubs; i++) {
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
  std::vector<int> d(1);
  cubatures_[0]->getAccuracy(d); degree_[0] = d[0];
  cubatures_[1]->getAccuracy(d); degree_[1] = d[0];

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
  std::vector<int> d(1);
  cubatures_[0]->getAccuracy(d); degree_[0] = d[0];
  cubatures_[1]->getAccuracy(d); degree_[1] = d[0];
  cubatures_[2]->getAccuracy(d); degree_[2] = d[0];

  dimension_ = cubatures_[0]->getDimension() + cubatures_[1]->getDimension() + cubatures_[2]->getDimension();
}



template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureTensor<Scalar,ArrayPoint,ArrayWeight>::CubatureTensor(Teuchos::RCP<CubatureDirect<Scalar,ArrayPoint,ArrayWeight> > cubature, int n) {
  cubatures_.resize(n);
  for (int i=0; i<n; i++) {
    cubatures_[i] = cubature;
  }

  std::vector<int> d(1);
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
  TEST_FOR_EXCEPTION( ( ( (int)cubPoints.size() < numCubPoints*cubDim ) || ( (int)cubWeights.size() < numCubPoints ) ),
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
