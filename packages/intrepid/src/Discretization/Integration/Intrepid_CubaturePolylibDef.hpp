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

/** \file   Intrepid_CubaturePolylibDef.hpp
    \brief  Definition file for the Intrepid::CubaturePolylib class.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

template <class Scalar, class ArrayPoint, class ArrayWeight>
CubaturePolylib<Scalar,ArrayPoint,ArrayWeight>::CubaturePolylib(int degree, EIntrepidPLPoly poly_type, Scalar alpha, Scalar beta) {
  TEST_FOR_EXCEPTION((degree < 0),
                     std::out_of_range,
                     ">>> ERROR (CubaturePolylib): No cubature rule implemented for the desired polynomial degree.");
  degree_    = degree;
  dimension_ = 1;
  poly_type_ = poly_type;
  alpha_     = alpha;
  beta_      = beta;
} // end constructor



template <class Scalar, class ArrayPoint, class ArrayWeight>
const char* CubaturePolylib<Scalar,ArrayPoint,ArrayWeight>::getName() const {
  return cubature_name_;
} // end getName



template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubaturePolylib<Scalar,ArrayPoint,ArrayWeight>::getDimension() const {
  return dimension_;
} // end dimension



template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubaturePolylib<Scalar,ArrayPoint,ArrayWeight>::getNumPoints() const {
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



template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubaturePolylib<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(std::vector<int> & accuracy) const {
  accuracy.assign(1, degree_);
} // end getAccuracy



template <class Scalar, class ArrayPoint, class ArrayWeight>
const char* CubaturePolylib<Scalar,ArrayPoint,ArrayWeight>::cubature_name_ = "INTREPID_CUBATURE_POLYLIB";



template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubaturePolylib<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint & cubPoints, ArrayWeight & cubWeights) const {
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
