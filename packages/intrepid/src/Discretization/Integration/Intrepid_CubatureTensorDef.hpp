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

template <class Scalar, class ArrayType>
CubatureTensor<Scalar,ArrayType>::CubatureTensor(std::vector< Teuchos::RCP<Cubature<Scalar,ArrayType> > > cubatures) {
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



template <class Scalar, class ArrayType>
CubatureTensor<Scalar,ArrayType>::CubatureTensor(Teuchos::RCP<CubatureDirect<Scalar,ArrayType> > cubature1,
                                                 Teuchos::RCP<CubatureDirect<Scalar,ArrayType> > cubature2) {
  cubatures_.resize(2);
  cubatures_[0] = cubature1;
  cubatures_[1] = cubature2;

  degree_.assign(2, 0);
  std::vector<int> d(1);
  cubatures_[0]->getAccuracy(d); degree_[0] = d[0];
  cubatures_[1]->getAccuracy(d); degree_[1] = d[0];

  dimension_ = cubatures_[0]->getDimension() + cubatures_[1]->getDimension();
}



template <class Scalar, class ArrayType>
CubatureTensor<Scalar,ArrayType>::CubatureTensor(Teuchos::RCP<CubatureDirect<Scalar,ArrayType> > cubature1,
                                                 Teuchos::RCP<CubatureDirect<Scalar,ArrayType> > cubature2,
                                                 Teuchos::RCP<CubatureDirect<Scalar,ArrayType> > cubature3) {
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



template <class Scalar, class ArrayType>
CubatureTensor<Scalar,ArrayType>::CubatureTensor(Teuchos::RCP<CubatureDirect<Scalar,ArrayType> > cubature, int n) {
  cubatures_.resize(n);
  for (int i=0; i<n; i++) {
    cubatures_[i] = cubature;
  }

  std::vector<int> d(1);
  cubatures_[0]->getAccuracy(d);
  degree_.assign(n,d[0]);

  dimension_ = cubatures_[0]->getDimension()*n;
}



template <class Scalar, class ArrayType>
void CubatureTensor<Scalar,ArrayType>::getCubature(ArrayType & cubPoints,
                                                   ArrayType & cubWeights) const {
  int numCubPoints = getNumPoints();
  int cubDim       = getDimension();
  // check size of cubPoints and cubWeights
  TEST_FOR_EXCEPTION( ( ( (int)cubPoints.size() < numCubPoints*cubDim ) || ( (int)cubWeights.size() < numCubPoints ) ),
                      std::out_of_range,
                      ">>> ERROR (CubatureTensor): Insufficient space allocated for cubature points or weights.");

  // set all weights to 1.0
  for (int i=0; i<numCubPoints; i++) {
      cubWeights(i) = (Scalar)1.0;
  }

  unsigned numCubs   = cubatures_.size();
  int globDimCounter = 0;
  int shift          = 1;
  for (unsigned i=0; i<numCubs; i++) {
    int numLocPoints = cubatures_[i]->getNumPoints();
    int locDim       = cubatures_[i]->getDimension();
    FieldContainer<Scalar> points(numLocPoints, locDim);
    FieldContainer<Scalar> weights(numLocPoints);
    cubatures_[i]->getCubature(points, weights);

    for (int j=0; j<numCubPoints; j++) {
      /* int itmp = ((j*shift) % numCubPoints) + (j / (numCubPoints/shift)); // equivalent, but numerically unstable */
      int itmp = (j % (numCubPoints/shift))*shift + (j / (numCubPoints/shift));
      for (int k=0; k<locDim; k++) {
        cubPoints(itmp , globDimCounter+k) = points(j % numLocPoints, k);
      }
      cubWeights( itmp ) *= weights(j % numLocPoints);
    }
    
    shift *= numLocPoints;
    globDimCounter += locDim;
  }
} // end getCubature



template <class Scalar, class ArrayType>
int CubatureTensor<Scalar,ArrayType>::getNumPoints() const {
  unsigned numCubs = cubatures_.size();
  int numCubPoints = 1;
  for (unsigned i=0; i<numCubs; i++) {
    numCubPoints *= cubatures_[i]->getNumPoints();
  }
  return numCubPoints;
} // end getNumPoints


template <class Scalar, class ArrayType>
int CubatureTensor<Scalar,ArrayType>::getDimension() const {
  return dimension_;
} // end dimension



template <class Scalar, class ArrayType>
void CubatureTensor<Scalar,ArrayType>::getAccuracy(std::vector<int> & degree) const {
  degree = degree_;
} // end getAccuracy

} // end namespace Intrepid
