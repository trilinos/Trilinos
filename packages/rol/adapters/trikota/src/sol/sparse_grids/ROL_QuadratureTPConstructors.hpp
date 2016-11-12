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

#ifndef ROL_QUADRATURETPCONSTRUCTORS_HPP
#define ROL_QUADRATURETPCONSTRUCTORS_HPP

namespace ROL {

// Build tensor product quadrature using default growth rule.
template <class Real> 
Quadrature<Real>::Quadrature(const int dimension,
                             const std::vector<int> &numPoints1D,
                             const std::vector<EQuadrature> &rule1D, 
                             const bool isNormalized) 
  : dimension_(dimension) {
  TEUCHOS_TEST_FOR_EXCEPTION((dimension!=(int)numPoints1D.size()||
		      dimension!=(int)rule1D.size()),std::out_of_range,
           ">>> ERROR (ROL::Quadrature): Dimension mismatch for inputs.");
  accuracy_.clear();
  std::vector<int> degree(1,0);
  Quadrature<Real> newRule(0,1);
  for (int i=0; i<dimension; i++) {
    // Compute 1D rules
    Quadrature<Real> rule1(rule1D[i],numPoints1D[i],isNormalized);
    rule1.getAccuracy(degree);
    accuracy_.push_back(degree[0]);
    // Build Tensor Rule
    newRule = kron_prod<Real>(newRule,rule1);
  }
  typename std::map<std::vector<Real>,int>::iterator it;
  int loc = 0;
  std::vector<Real> node(dimension,0.0);
  for (it=newRule.begin(); it!=newRule.end(); it++) {
    node = it->first;
    addPointAndWeight(node,newRule.getWeight(node),loc);
    loc++;
  }
  if (isNormalized) {
    normalize();
  }
}

// Build tensor product quadrature using user-input growth rule.
template <class Real> 
Quadrature<Real>::Quadrature(const int dimension,
                             const std::vector<int> &numPoints1D,
                             const std::vector<EQuadrature> &rule1D, 
                             const std::vector<EGrowth> &growth1D,
                             const bool isNormalized) 
  : dimension_(dimension) {
  TEUCHOS_TEST_FOR_EXCEPTION((dimension!=(int)numPoints1D.size()||
		      dimension!=(int)rule1D.size()||
		      dimension!=(int)growth1D.size()),std::out_of_range,
           ">>> ERROR (Quadrature): Dimension mismatch for inputs.");
  accuracy_.clear();
  accuracy_.resize(dimension);
  std::vector<int> degree(1);
  Quadrature<Real> newRule(0,1);
  for (int i=0; i<dimension; i++) {
    // Compute 1D rules
    int numPoints = growthRule1D(numPoints1D[i],growth1D[i],rule1D[i]);
    Quadrature<Real> rule1(rule1D[i],numPoints,isNormalized);
    rule1.getAccuracy(degree);
    accuracy_.push_back(degree[0]);
    // Build Tensor Rule
    newRule = kron_prod<Real>(newRule,rule1);
  }
  typename std::map<std::vector<Real>,int>::iterator it;
  int loc = 0;
  std::vector<Real> node;
  for (it=newRule.begin(); it!=newRule.end(); it++) {
    node = it->first;
    addPointAndWeight(node,newRule.getWeight(node),loc);
    loc++;
  }
  if (isNormalized) {
    normalize();
  }
}

// Build (adaptive) sparse-grid quadrature rule.
template <class Real> 
Quadrature<Real>::Quadrature(const int dimension,
                             const int maxNumPoints,
                             const std::vector<EQuadrature> &rule1D, 
                             const std::vector<EGrowth> &growth1D,
                             const bool isNormalized,
                             const bool adaptive)
  : dimension_(dimension) {
  TEUCHOS_TEST_FOR_EXCEPTION((dimension!=(int)rule1D.size()||
		      dimension!=(int)growth1D.size()),std::out_of_range,
            ">>> ERROR (Quadrature): Dimension mismatch for inputs.");
  if ( adaptive ) {
    buildInitial(dimension,maxNumPoints,rule1D,growth1D,isNormalized);
  }
  else {
    buildSGMGA(dimension,maxNumPoints,rule1D,growth1D);
  }
  if (isNormalized) {
    normalize();
  }
} 

// Build (adaptive) sparse-grid quadrature rule.
template <class Real> 
Quadrature<Real>::Quadrature(const QuadratureInfo &info) {
  dimension_ = info.dim;
  int maxNumPoints = info.maxLevel-1;
  std::vector<EQuadrature> rule1D = info.rule1D;
  std::vector<EGrowth> growth1D = info.growth1D;
  bool isNormalized = info.normalized;
  bool adaptive = info.adaptive;  
  TEUCHOS_TEST_FOR_EXCEPTION(((dimension_!=(int)rule1D.size()) ||
		              (dimension_!=(int)growth1D.size())),
                             std::out_of_range,
    ">>> ERROR (Quadrature): Dimension mismatch for inputs.");
  if ( adaptive ) {
    buildInitial(dimension_,maxNumPoints,rule1D,growth1D,isNormalized);
  }
  else {
    buildSGMGA(dimension_,maxNumPoints,rule1D,growth1D);
  }
  if (isNormalized) {
    normalize();
  }
} 

// Build quadrature rule from input files.
template<class Real>
Quadrature<Real>::Quadrature(const char* SGinfo,
                             const char* SGdata,
                             const bool isNormalized) {
  // Open SGdata
  std::fstream info; info.open(SGinfo);
  Real buf;

  // Get Dimension Number
  info >> buf;
  int dim_num = (int)buf;
  dimension_ = dim_num;

  // Get Number of Cubature Points/Weights
  info >> buf;
  int point_num = (int)buf;

  // Close SGinfo
  info.close();

  // Open SGdata File
  std::fstream file; file.open(SGdata);

  // Get Cubature Points/Weights
  Real weight = 0.0;
  std::vector<Real> node(dim_num);
  int loc = 0;
  for ( int i = 0; i < point_num; i++ ) {
    // Get Cubature Point
    for ( int j = 0; j < dim_num; j++ ) {
      file >> buf;
      node[j] = buf;
    }
    // Get Cubature Weight
    file >> buf;
    weight = buf;
    // Insert Point/Weight into Sparse Grid Storage
    addPointAndWeight(node,weight,loc);
    loc++;
    // Move to Next Line in SGdata
    file.ignore(256,'\n');
  }
  file.close();
  if (isNormalized) {
    normalize();
  }
}
} // end ROL namespace

#endif
