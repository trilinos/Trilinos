// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file   ROL_QuadratureDef.hpp
    \brief  Definition file for the ROL::Quadrature class.
    \author Created by D. Kouri and D. Ridzal.
*/

namespace ROL {

template <class Real> 
Quadrature<Real>::Quadrature( int degree, EROLBurkardt rule, bool isNormalized ) 
  : dimension_(1) {
  TEUCHOS_TEST_FOR_EXCEPTION((degree < 0),std::out_of_range,
    ">>> ERROR (Quadrature): No rule implemented for desired polynomial degree.");
  accuracy_.clear(); 
  accuracy_.push_back(degree);

  int numPoints = (degree+1)/2+1;
  if (rule==BURK_CLENSHAWCURTIS||rule==BURK_FEJER2) {
    numPoints = degree+1;
  }
  else if (rule==BURK_TRAPEZOIDAL) {
    numPoints = 2;
  }
  else if (rule==BURK_PATTERSON) {
    int l = 0, o = (degree-0.5)/1.5;
    for (int i=0; i<8; i++) {
      l = (int)pow(2.0,(double)i+1.0)-1;
      if (l>=o) {
	numPoints = l;
	break;
      }
    }
  }
  else if (rule==BURK_GENZKEISTER) {
    int o_ghk[8] = {1,3,9,19,35,37,41,43}; 
    int o = (degree-0.5)/1.5;
    for (int i=0; i<8; i++) {
      if (o_ghk[i]>=o) {
	numPoints = o_ghk[i];
	break;
      }
    }
  }

  std::vector<Real> points(numPoints), weights(numPoints);
  switch(rule) {
    case BURK_CHEBYSHEV1: // Gauss-Chebyshev Type 1
      ROLBurkardtRules::chebyshev1_compute<Real>(numPoints,&points[0],&weights[0]); break;
    case BURK_CHEBYSHEV2: // Gauss-Chebyshev Type 2
      ROLBurkardtRules::chebyshev2_compute<Real>(numPoints,&points[0],&weights[0]); break;
    case BURK_CLENSHAWCURTIS: // Clenshaw-Curtis    
      ROLBurkardtRules::clenshaw_curtis_compute<Real>(numPoints,&points[0],&weights[0]); break;
    case BURK_FEJER2: // Fejer Type 2
      ROLBurkardtRules::fejer2_compute<Real>(numPoints,&points[0],&weights[0]); break;
    case BURK_LEGENDRE: // Gauss-Legendre
      ROLBurkardtRules::legendre_compute<Real>(numPoints,&points[0],&weights[0]); break;
    case BURK_PATTERSON: // Gauss-Patterson
      ROLBurkardtRules::patterson_lookup<Real>(numPoints,&points[0],&weights[0]); break;
    case BURK_TRAPEZOIDAL: // Trapezoidal Rule
      ROLBurkardtRules::trapezoidal_compute<Real>(numPoints,&points[0],&weights[0]); break;
    case BURK_HERMITE: // Gauss-Hermite
      ROLBurkardtRules::hermite_compute<Real>(numPoints,&points[0],&weights[0]); break;
    case BURK_GENZKEISTER: // Hermite-Genz-Keister
      ROLBurkardtRules::hermite_genz_keister_lookup<Real>(numPoints,&points[0],&weights[0]); break;
    case BURK_LAGUERRE: // Gauss-Laguerre
      ROLBurkardtRules::laguerre_compute<Real>(numPoints,&points[0],&weights[0]); break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::out_of_range,
        ">>> ERROR (Quadrature): No such 1D quadrature rule."); break;
  }

  if (isNormalized) {
    Real sum = 0.0;
    for (int i=0; i<numPoints; i++) {
      sum += weights[i];
    }
    for (int i=0; i<numPoints; i++) {
      weights[i] /= sum;
    }
  }
  else {
    if (rule==BURK_GENZKEISTER && numPoints >= 37) {
      for (int i=0; i<numPoints; i++) {
	weights[i] *= sqrt(M_PI);
      }
    }
  }

  std::vector<Real> point(1,0.0);
  for (int i=0; i<numPoints; i++) {
    point[0] = points[i];
    addPointAndWeight(point,weights[i],i);
  }
} // end constructor
  
template <class Real> 
Quadrature<Real>::Quadrature( EROLBurkardt rule, int numPoints, bool isNormalized )
  : dimension_(1) {
  TEUCHOS_TEST_FOR_EXCEPTION((numPoints < 0),std::out_of_range, 
     ">>> ERROR (Quadrature): No rule implemented for desired number of points.");

  accuracy_.clear();
  std::vector<Real> points(numPoints), weights(numPoints);
  if (rule==BURK_CHEBYSHEV1) { // Gauss-Chebyshev Type 1
    accuracy_.push_back(2*numPoints-1);
    ROLBurkardtRules::chebyshev1_compute<Real>(numPoints,&points[0],&weights[0]);
  }
  else if (rule==BURK_CHEBYSHEV2) { // Gauss-Chebyshev Type 2
    accuracy_.push_back(2*numPoints-1);
    ROLBurkardtRules::chebyshev2_compute<Real>(numPoints,&points[0],&weights[0]);
  }
  else if (rule==BURK_CLENSHAWCURTIS) { // Clenshaw-Curtis    
    accuracy_.push_back(numPoints-1);
    ROLBurkardtRules::clenshaw_curtis_compute<Real>(numPoints,&points[0],&weights[0]);
  }
  else if (rule==BURK_FEJER2) { // Fejer Type 2
    accuracy_.push_back(numPoints-1);
    ROLBurkardtRules::fejer2_compute<Real>(numPoints,&points[0],&weights[0]);
  }
  else if (rule==BURK_LEGENDRE) { // Gauss-Legendre
    accuracy_.push_back(2*numPoints-1);
    ROLBurkardtRules::legendre_compute<Real>(numPoints,&points[0],&weights[0]);
  }
  else if (rule==BURK_PATTERSON) { // Gauss-Patterson
    bool correctNumPoints = false;
    for (int i=0; i<8; i++) {
      int l = (int)pow(2.0,(double)i+1.0)-1;
      if (numPoints==l) {
	correctNumPoints = true;
	break;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION((correctNumPoints==false),std::out_of_range,
	">>> ERROR (Quadrature): Number of points must be numPoints = 1, 3, 7, 15, 31, 63, 127, 255.");
    Real degree = 1.5*(double)numPoints+0.5;
    accuracy_.push_back((int)degree);
    ROLBurkardtRules::patterson_lookup<Real>(numPoints,&points[0],&weights[0]);
  }
  else if (rule==BURK_TRAPEZOIDAL) { // Trapezoidal Rule
    accuracy_.push_back(2);
    ROLBurkardtRules::trapezoidal_compute<Real>(numPoints,&points[0],&weights[0]);
  }
  else if (rule==BURK_HERMITE) { // Gauss-Hermite
    accuracy_.push_back(2*numPoints-1);
    ROLBurkardtRules::hermite_compute<Real>(numPoints,&points[0],&weights[0]);
  }
  else if (rule==BURK_GENZKEISTER) { // Hermite-Genz-Keister
    bool correctNumPoints = false;
    int o_ghk[8] = {1,3,9,19,35,37,41,43};
    for (int i=0; i<8; i++) {
      if (o_ghk[i]==numPoints) {
	correctNumPoints = true;
	break;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION((correctNumPoints==false),std::out_of_range,
       ">>> ERROR (Quadrature): Number of points must be numPoints = 1, 3, 9, 35, 37, 41, 43.");
    Real degree = 1.5*(double)numPoints+0.5;
    accuracy_.push_back((int)degree);
    ROLBurkardtRules::hermite_genz_keister_lookup<Real>(numPoints,&points[0],&weights[0]);
  }
  else if (rule==BURK_LAGUERRE) { // Gauss-Laguerre
    accuracy_.push_back(2*numPoints-1);
    ROLBurkardtRules::laguerre_compute<Real>(numPoints,&points[0],&weights[0]);
  }
  
  if (isNormalized) {
    Real sum = 0.0;
    for (int i=0; i<numPoints; i++) {
      sum += weights[i];
    }
    for (int i=0; i<numPoints; i++) {
      weights[i] /= sum;
    }
  }
  else {
    if (rule==BURK_GENZKEISTER && numPoints >= 37) {
      for (int i=0; i<numPoints; i++) {
	weights[i] *= sqrt(M_PI);
      }
    }
  }

  std::vector<Real> point(1,0.0);
  for (int i=0; i<numPoints; i++) {
    point[0] = points[i];
    addPointAndWeight(point,weights[i],i);
  }
} // end constructor
  
template <class Real>
Quadrature<Real>::Quadrature( std::vector<Real>& points, std::vector<Real>& weights) 
  : dimension_(1) {
  int size = (int)weights.size();
  TEUCHOS_TEST_FOR_EXCEPTION(((int)points.size()!=size),std::out_of_range,
	     ">>> ERROR (Quadrature): Input dimension mismatch.");
  std::vector<Real> point(1,0.0);
  for (int i=0; i<size; i++) {
    point[0] = points[i];
    addPointAndWeight(point,weights[i],i);
  }
}

} // ROL namespace
