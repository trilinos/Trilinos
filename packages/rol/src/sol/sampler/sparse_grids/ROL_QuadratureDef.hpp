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
Quadrature<Real>::Quadrature(int dimension) : dimension_(dimension) {
  points_.clear(); weights_.clear(); accuracy_.clear();
  accuracy_.resize(dimension);
}

template <class Real>
void Quadrature<Real>::addPointAndWeight(std::vector<Real> point, Real weight, int loc) {
  this->points_.insert(std::pair<std::vector<Real>,int>(point,loc));
  this->weights_.push_back(weight);
}

template <class Real>
int Quadrature<Real>::getNumPoints() const {
  return this->weights_.size();
} 

template <class Real>
void Quadrature<Real>::getAccuracy(std::vector<int> & accuracy) const {
  accuracy = this->accuracy_;
}

template <class Real>
void Quadrature<Real>::getCubature(std::vector<std::vector<Real> >& points, 
                                   std::vector<Real>& weights) const {
  points.resize(this->weights_.size());
  weights.resize(this->weights_.size());
  typename std::map<std::vector<Real>,int>::const_iterator it;
  for (it=(this->points_).begin(); it!=(this->points_).end();it++) {
    points[it->second] = it->first;
    weights[it->second] = (this->weights_)[it->second];
  }
}

template <class Real>
int Quadrature<Real>::getDimension() const {
  return this->dimension_;
}

template <class Real> 
typename std::map<std::vector<Real>,int>::iterator Quadrature<Real>::begin() {
  return (this->points_).begin();
}

template <class Real> 
typename std::map<std::vector<Real>,int>::iterator Quadrature<Real>::end() {
  return (this->points_).end();
}

template <class Real> 
void Quadrature<Real>::insert(typename std::map<std::vector<Real>,int>::iterator it, 
                              std::vector<Real> point, Real weight) {
  (this->points_).insert(it,std::pair<std::vector<Real>,int>(point,(int)(this->points_).size()));
  (this->weights_).push_back(weight);
}

template <class Real> 
std::vector<Real> Quadrature<Real>::getNode(typename std::map<std::vector<Real>,int>::iterator it) {
  return it->first;
}

template <class Real> 
Real Quadrature<Real>::getWeight(int node) { 
  return (this->weights_)[node];
}

template <class Real> 
Real Quadrature<Real>::getWeight(std::vector<Real> point) {
  return (this->weights_)[(this->points_)[point]];
}

template <class Real>
void Quadrature<Real>::update(Real alpha, Quadrature<Real> &rule) {
  if ( rule.begin() != rule.end() ) {  
    // Initialize an iterator on std::map<std::vector<Real>,Real>
    typename std::map<std::vector<Real>,int>::iterator it;
    typename std::map<std::vector<Real>,int>::iterator it2;
  
    // Get location of last point and weight 
    int loc = this->points_.size();
 
    // Intersection of rule1 and rule2 and set difference rule2 \ rule1
    for ( it2=rule.begin(); it2!=rule.end(); it2++ ) {
      it = (this->points_).find(it2->first);
      if ( it != (this->points_).end() ) {
        (this->weights_)[it->second] += alpha*rule.getWeight(it2->second);
      } 
      else {
        (this->points_).insert(std::pair<std::vector<Real>,int>(it2->first,loc));
        (this->weights_).push_back(alpha*rule.getWeight(it2->second));
        loc++;
      }
    }
  }
}

template <class Real>
void Quadrature<Real>::normalize(){
  Real sum = 0.0;
  typename std::vector<Real>::iterator it;
  for (it=weights_.begin(); it!=weights_.end(); it++) {
    sum += *it;
  }
  for (it=weights_.begin(); it!=weights_.end(); it++) {
    *it /= sum;
  }
}

} // end ROL namespace
