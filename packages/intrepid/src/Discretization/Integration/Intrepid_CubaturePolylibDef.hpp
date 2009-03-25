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

/** \file   Intrepid_CubaturePolylibDef.hpp
    \brief  Definition file for the Intrepid::CubaturePolylib class.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

template <class Scalar, class ArrayType>
CubaturePolylib<Scalar,ArrayType>::CubaturePolylib(int degree, EIntrepidPLPoly poly_type, Scalar alpha, Scalar beta) {
  TEST_FOR_EXCEPTION((degree < 0),
                     std::out_of_range,
                     ">>> ERROR (CubaturePolylib): No cubature rule implemented for the desired polynomial degree.");
  degree_    = degree;
  dimension_ = 1;
  poly_type_ = poly_type;
  alpha_     = alpha;
  beta_      = beta;
} // end constructor



template <class Scalar, class ArrayType>
const char* CubaturePolylib<Scalar,ArrayType>::getName() const {
  return cubature_name_;
} // end getName



template <class Scalar, class ArrayType>
int CubaturePolylib<Scalar,ArrayType>::getDimension() const {
  return dimension_;
} // end dimension



template <class Scalar, class ArrayType>
int CubaturePolylib<Scalar,ArrayType>::getNumPoints() const {
  int np = 0;
  switch (poly_type_) {
    case PL_GAUSS:
      np = (degree_+(int)2)/(int)2;
      break;
    case PL_GAUSS_RADAU_LEFT:
    case PL_GAUSS_RADAU_RIGHT:
      if (degree_ == 0)
        np = 2;
      else
        np = (degree_+(int)3)/(int)2;
      break;
    case PL_GAUSS_LOBATTO:
      np = (degree_+(int)4)/(int)2;
      break;
    default:
      TEST_FOR_EXCEPTION((1),
                         std::invalid_argument,
                         ">>> ERROR (CubaturePolylib): Unknown point type argument.");
  }
  return np;
} // end getNumPoints



template <class Scalar, class ArrayType>
void CubaturePolylib<Scalar,ArrayType>::getAccuracy(std::vector<int> & accuracy) const {
  accuracy.assign(1, degree_);
} // end getAccuracy



template <class Scalar, class ArrayType>
const char* CubaturePolylib<Scalar,ArrayType>::cubature_name_ = "INTREPID_CUBATURE_POLYLIB";



template <class Scalar, class ArrayType>
void CubaturePolylib<Scalar,ArrayType>::getCubature(ArrayType & cubPoints, ArrayType & cubWeights) const {
  int numCubPoints = getNumPoints();
  int cellDim      = getDimension();
  // check size of cubPoints and cubWeights
  TEST_FOR_EXCEPTION( ( ( (int)cubPoints.size() < numCubPoints*cellDim ) || ( (int)cubWeights.size() < numCubPoints ) ),
                      std::out_of_range,
                      ">>> ERROR (CubatureDirect): Insufficient space allocated for cubature points or weights.");

  // temporary storage
  FieldContainer<Scalar> z(numCubPoints);
  FieldContainer<Scalar> w(numCubPoints);

  // run Polylib routines
  switch (poly_type_) {
    case PL_GAUSS:
      IntrepidPolylib::zwgj(&z[0], &w[0], numCubPoints, alpha_, beta_);
      break;
    case PL_GAUSS_RADAU_LEFT:
      IntrepidPolylib::zwgrjm(&z[0], &w[0], numCubPoints, alpha_, beta_);
      break;
    case PL_GAUSS_RADAU_RIGHT:
      IntrepidPolylib::zwgrjp(&z[0], &w[0], numCubPoints, alpha_, beta_);
      break;
    case PL_GAUSS_LOBATTO:
      IntrepidPolylib::zwglj(&z[0], &w[0], numCubPoints, alpha_, beta_);
      break;
    default:
      TEST_FOR_EXCEPTION((1),
                         std::invalid_argument,
                         ">>> ERROR (CubaturePolylib): Unknown point type argument.");
  }

  // fill input arrays
  for (int pointId = 0; pointId < numCubPoints; pointId++) {
    for (int dim = 0; dim < cellDim; dim++) {
      cubPoints(pointId,dim) = z[pointId];
    }
    cubWeights(pointId) = w[pointId];
  }
} // end getCubature


} // end namespace Intrepid
