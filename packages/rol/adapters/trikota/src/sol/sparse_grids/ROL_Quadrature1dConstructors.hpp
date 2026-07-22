// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_QUADRATURE1DCONSTRUCTORS_HPP
#define ROL_QUADRATURE1DCONSTRUCTORS_HPP

/** \file   ROL_Quadrature1DConstructors.hpp
    \brief  One-dimensional quadrature constructors for the ROL::Quadrature class.
    \author Created by D. Kouri and D. Ridzal.
*/

#include "TriKota_ConfigDefs.hpp"
#include "sandia_rules.hpp"

namespace ROL {

template <class Real> 
Quadrature<Real>::Quadrature(const int degree,
                             const EQuadrature rule,
                             const bool isNormalized) 
  : dimension_(1) {
  ROL_TEST_FOR_EXCEPTION((degree < 0),std::out_of_range,
    ">>> ERROR (Quadrature): No rule implemented for desired polynomial degree.");
  accuracy_.clear(); 
  accuracy_.push_back(degree);

  int numPoints = (degree+1)/2+1;
  if (rule==QUAD_CLENSHAWCURTIS||rule==QUAD_FEJER2) {
    numPoints = degree+1;
  }
  else if (rule==QUAD_TRAPEZOIDAL) {
    numPoints = 2;
  }
  else if (rule==QUAD_PATTERSON) {
    int l = 0, o = (degree-0.5)/1.5;
    for (int i=0; i<8; i++) {
      l = (int)pow(2.0,(double)i+1.0)-1;
      if (l>=o) {
	numPoints = l;
	break;
      }
    }
  }
  else if (rule==QUAD_GENZKEISTER) {
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
    case QUAD_CHEBYSHEV1: // Gauss-Chebyshev Type 1
      webbur::chebyshev1_compute(numPoints,&points[0],&weights[0]); break;
    case QUAD_CHEBYSHEV2: // Gauss-Chebyshev Type 2
      webbur::chebyshev2_compute(numPoints,&points[0],&weights[0]); break;
    case QUAD_CLENSHAWCURTIS: // Clenshaw-Curtis    
      webbur::clenshaw_curtis_compute(numPoints,&points[0],&weights[0]); break;
    case QUAD_FEJER2: // Fejer Type 2
      webbur::fejer2_compute(numPoints,&points[0],&weights[0]); break;
    case QUAD_LEGENDRE: // Gauss-Legendre
      webbur::legendre_compute(numPoints,&points[0],&weights[0]); break;
    case QUAD_PATTERSON: // Gauss-Patterson
      webbur::patterson_lookup(numPoints,&points[0],&weights[0]); break;
    case QUAD_TRAPEZOIDAL: { // Trapezoidal Rule
      Real zero(0), one(1), two(2);
      if (numPoints==1) {
        points[0] = zero;
        weights[0] = two;
      }
      else {
        Real h = one/(static_cast<Real>(numPoints)-one);
        for (int i = 0; i < numPoints; ++i) {
          points[i] = -one + static_cast<Real>(i)*two*h;
          if (i==0 || i==numPoints-1) {
            weights[i] = h;
          }
          else {
            weights[i] = two*h;
          }
        }
      }
      break;
    }
    case QUAD_HERMITE: // Gauss-Hermite
      webbur::hermite_compute(numPoints,&points[0],&weights[0]); break;
    case QUAD_GENZKEISTER: // Hermite-Genz-Keister
      webbur::hermite_genz_keister_lookup(numPoints,&points[0],&weights[0]); break;
    case QUAD_LAGUERRE: // Gauss-Laguerre
      webbur::laguerre_compute(numPoints,&points[0],&weights[0]); break;
    default:
      ROL_TEST_FOR_EXCEPTION(true,std::out_of_range,
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
    if (rule==QUAD_GENZKEISTER && numPoints >= 37) {
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
Quadrature<Real>::Quadrature(const EQuadrature rule,
                             const int numPoints,
                             const bool isNormalized)
  : dimension_(1) {
  ROL_TEST_FOR_EXCEPTION((numPoints < 0),std::out_of_range, 
     ">>> ERROR (Quadrature): No rule implemented for desired number of points.");

  accuracy_.clear();
  std::vector<Real> points(numPoints), weights(numPoints);
  if (rule==QUAD_CHEBYSHEV1) { // Gauss-Chebyshev Type 1
    accuracy_.push_back(2*numPoints-1);
    webbur::chebyshev1_compute(numPoints,&points[0],&weights[0]);
  }
  else if (rule==QUAD_CHEBYSHEV2) { // Gauss-Chebyshev Type 2
    accuracy_.push_back(2*numPoints-1);
    webbur::chebyshev2_compute(numPoints,&points[0],&weights[0]);
  }
  else if (rule==QUAD_CLENSHAWCURTIS) { // Clenshaw-Curtis    
    accuracy_.push_back(numPoints-1);
    webbur::clenshaw_curtis_compute(numPoints,&points[0],&weights[0]);
  }
  else if (rule==QUAD_FEJER2) { // Fejer Type 2
    accuracy_.push_back(numPoints-1);
    webbur::fejer2_compute(numPoints,&points[0],&weights[0]);
  }
  else if (rule==QUAD_LEGENDRE) { // Gauss-Legendre
    accuracy_.push_back(2*numPoints-1);
    webbur::legendre_compute(numPoints,&points[0],&weights[0]);
  }
  else if (rule==QUAD_PATTERSON) { // Gauss-Patterson
    bool correctNumPoints = false;
    for (int i=0; i<8; i++) {
      int l = (int)pow(2.0,(double)i+1.0)-1;
      if (numPoints==l) {
	correctNumPoints = true;
	break;
      }
    }
    ROL_TEST_FOR_EXCEPTION((correctNumPoints==false),std::out_of_range,
	">>> ERROR (Quadrature): Number of points must be numPoints = 1, 3, 7, 15, 31, 63, 127, 255.");
    Real degree = 1.5*(double)numPoints+0.5;
    accuracy_.push_back((int)degree);
    webbur::patterson_lookup(numPoints,&points[0],&weights[0]);
  }
  else if (rule==QUAD_TRAPEZOIDAL) { // Trapezoidal Rule
    accuracy_.push_back(2);
    Real zero(0), one(1), two(2);
    if (numPoints==1) {
      points[0] = zero;
      weights[0] = two;
    }
    else {
      Real h = one/(static_cast<Real>(numPoints)-one);
      for (int i = 0; i < numPoints; ++i) {
        points[i] = -one + static_cast<Real>(i)*two*h;
        if (i==0 || i==numPoints-1) {
          weights[i] = h;
        }
        else {
          weights[i] = two*h;
        }
      }
    }
  }
  else if (rule==QUAD_HERMITE) { // Gauss-Hermite
    accuracy_.push_back(2*numPoints-1);
    webbur::hermite_compute(numPoints,&points[0],&weights[0]);
  }
  else if (rule==QUAD_GENZKEISTER) { // Hermite-Genz-Keister
    bool correctNumPoints = false;
    int o_ghk[8] = {1,3,9,19,35,37,41,43};
    for (int i=0; i<8; i++) {
      if (o_ghk[i]==numPoints) {
	correctNumPoints = true;
	break;
      }
    }
    ROL_TEST_FOR_EXCEPTION((correctNumPoints==false),std::out_of_range,
       ">>> ERROR (Quadrature): Number of points must be numPoints = 1, 3, 9, 35, 37, 41, 43.");
    Real degree = 1.5*(double)numPoints+0.5;
    accuracy_.push_back((int)degree);
    webbur::hermite_genz_keister_lookup(numPoints,&points[0],&weights[0]);
  }
  else if (rule==QUAD_LAGUERRE) { // Gauss-Laguerre
    accuracy_.push_back(2*numPoints-1);
    webbur::laguerre_compute(numPoints,&points[0],&weights[0]);
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
    if (rule==QUAD_GENZKEISTER && numPoints >= 37) {
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
Quadrature<Real>::Quadrature(const std::vector<Real> &points,
                             const std::vector<Real> &weights) 
  : dimension_(1) {
  int size = (int)weights.size();
  ROL_TEST_FOR_EXCEPTION(((int)points.size()!=size),std::out_of_range,
	     ">>> ERROR (Quadrature): Input dimension mismatch.");
  std::vector<Real> point(1);
  for (int i=0; i<size; i++) {
    point[0] = points[i];
    addPointAndWeight(point,weights[i],i);
  }
}

} // ROL namespace

#endif
