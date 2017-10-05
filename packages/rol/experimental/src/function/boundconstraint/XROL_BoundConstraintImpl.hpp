
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

#pragma once

#include "XROL_VectorTraits.hpp"


namespace XROL {


template<class V>
BoundConstraint<V>::BoundConstraint( std::unique_ptr<V> x_lo, 
                                     std::unique_ptr<V> x_up,
                                     magnitude_t<V> scale ) :
  x_lo_(std::move(x_lo)), x_up_(std::move(x_up)), 
  mask_(std::move(clone(x_lo)),  
  scale_(scale), Lactivated_(true), Uactivated_(true) {

  // Create a minimum function for this element type
  auto min = make_min(*x_lo_);
  auto diff = []( auto u, auto l ) { return u-l; };

  // Compute half the minimum distance between the upper and lower bounds
  min_diff_ = Elementwise::eval_function_and_reduce(min,diff,*x_up_,*x_lo_);
  min_diff_ *= 0.5;

}

template<class V> 
BoundConstraint<V>::BoundConstraint( const V& x, 
                                     bool isLower, 
                                     magnitude_t<V> scale ) :
  x_lo_(std::move(clone(x))), 
  x_up_(std::move(clone(x))),
  scale_(scale),
  Lactivated_( isLower),
  Uactivated_(!isLower) {

  if( isLower ) {
    set(*x_lo_,x);
    fill(*x_up_,INF_); // Set all elements to "infinity"
  else {
    set(*x_up_,x);
    fill(*x_lo_,NINF_); // Set all elements to "negative infinity"
  }
}

template<class V>
BoundConstraint<V>::~BoundConstraint() {}


template<class V>
void BoundConstraint<V>::project( V& x ) const {
  // Scalar clipping function
  auto clip = []( auto ue, auto xe, auto le ) {
    return std::min( ue, std::max(xe,le) );
  };

  // x_i = min( u_i, max( l_i, x_i ) ) for all i
  Elementwise::eval_function(x,clip,*x_up_,x,*x_lo_);
}

/* Make vector strictly feasible

   x_i = min( umod_i, max( lmod_i, x_i ) ) for all i

   where  
*/
template<class V>
void BoundConstraint<V>::projectInterior( V& x ) const {

  auto eps = std::numeric_limits<magnitude_t<V>>::epsilon();
  auto tol = 100*eps;
  magnitude_t<V> one(1);
}

template<class V>
void BoundConstraint<V>::pruneUpperActive( V& v, 
                                           const V& x,  
                                           magnitude_t<V> eps ) {
  auto epsn(std::min(scale_*eps,min_diff_)); 
  auto active = [epsn]( auto xe, auto ye ) { return (y<=epsn) ? 0 : x; };

}

template<class V>
void BoundConstraint<V>::pruneUpperActive( V& v, 
                                           const dual_t<V> &g, 
                                           const V& x,  
                                           magnitude_t<V> eps ) {
  auto epsn(std::min(scale_*eps,min_diff_)); 

}

template<class V>
void BoundConstraint<V>::pruneLowerActive( V& v, 
                                           const V& x,  
                                           magnitude_t<V> eps ) {
  auto epsn(std::min(scale_*eps,min_diff_)); 

}

template<class V>
void BoundConstraint<V>::pruneLowerActive( V& v, 
                                           const dual_t<V> &g, 
                                           const V& x,  
                                           magnitude_t<V> eps ) {
  auto epsn(std::min(scale_*eps,min_diff_)); 

}

template<class V>
const V& BoundConstraint<V>::getLowerBound( void ) const {
  return *x_lo_;
}

template<class V>
const V& BoundConstraint<V>::getUpperBound( void ) const {
  return *x_up_;
}

bool isFeasible( const V& v ) const {
  auto is_between = []( auto l, auto v, auto u ) {
    return ( (l<v) && (v<u) ) ? magnitude_t<V>(1) : magnitude_t<V>(0);
  }; 
  auto all_true = make_product(v);
  return (Elementwise::eval_function_and_reduce(all_true,is_between,*x_lo_,v,*x_up_) != 0);
}

} // namespace XROL

