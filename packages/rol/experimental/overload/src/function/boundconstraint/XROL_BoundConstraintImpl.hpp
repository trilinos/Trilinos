
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
                                     Real scale ) :
  x_lo_(std::move(x_lo)), x_up_(std::move(x_up)), 
  mask_(std::move(clone(x_lo))),  
  scale_(scale), Lactivated_(true), Uactivated_(true) {

  // Create a minimum function for this element type
  auto min = make_min(*x_lo_);
  auto diff = []( auto u, auto l ) { return u-l; };

  // Compute half the minimum distance between the upper and lower bounds
  min_diff_ = Elementwise::eval_function_and_reduce(min,diff,*x_up_,*x_lo_);
  min_diff_ *= 0.5;

}

template<class V> 
BoundConstraint<V>::BoundConstraint( const V& x, bool isLower, Real scale ) :
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
    return std::min( le, std::max(xe,ue) );
  };

  /**  \f[ x_i = \min( u_i, \max( l_i, x_i ) ) \; \forall i \f] */
  Elementwise::eval_function(x,clip,*x_lo_,x,*x_up_);
}


template<class V>
void BoundConstraint<V>::projectInterior( V& x ) const {

  /* Same as project with modifications to lower and upper bounds: 
     \f[ x_i = \min( \tilde u_i, \max( \tilde l_i, x_i ) ) \; \forall i \f]
   */

  auto eps = std::numeric_limits<Real>::epsilon();
  auto tol = 100*eps;
  Real one(1);

  // Offset is positive for lower and negative for upper bounds
  // 
  auto mod = [eps,] ( auto offset, auto val ) { 
    return (val == 0) ? offset : (1+signbit(val)*offset)*val;
  };

  auto clip = [tol,&mod]( auto le, auto xe, auto ue ) {
    return std::min( mod(tol,le), std::max( xe, mod(-tol,ue) );
  };

  Elementwise::eval_function(x,clip,*x_lo_,x,*x_up_);

}



template<class V>
void BoundConstraint<V>::pruneUpperActive( V& v, const V& x, Real eps ) const {
  auto epsn{std::min(scale_*eps,min_diff_)}; 
  auto f = [epsn]( auto ve, auto xe, auto ue) { 
    return ( (ue-xe)<=epsn ) ? 0 : ve; 
  };
  Elementwise::eval_function(v,f,x,*x_up_);
}

template<class V>
void BoundConstraint<V>::pruneUpperActive( V& v, const dual_t<V>& g, const V& x, Real eps ) const {
  auto epsn{std::min(scale_*eps,min_diff_)}; 
  auto f = [epsn]( auto ve, auto xe, auto ue, auto ge ) { 
    return ( ((ue-xe)<=epsn) && ge<0 ) ? 0 : ve; 
  };
  Elementwise::eval_function(v,f,x,*x_up_,g);
}

template<class V>
void BoundConstraint<V>::pruneLowerActive( V& v, const V& x, Real eps ) const {
  auto epsn{std::min(scale_*eps,min_diff_)}; 
  auto f = [epsn]( auto ve, auto xe, auto le) { 
    return ( (xe-le)<=epsn ) ? 0 : ve; 
  };
  Elementwise::eval_function(v,f,x,*x_lo_);
}

template<class V>
void BoundConstraint<V>::pruneLowerActive( V& v, const dual_t<V>& g, const V& x, Real eps ) const {
  auto epsn{std::min(scale_*eps,min_diff_)}; 
  auto f = [epsn]( auto ve, auto xe, auto le, auto ge ) { 
    return ( ((xe-le)<=epsn) && ge>0 ) ? 0 : ve; 
  };
  Elementwise::eval_function(v,f,x,*x_lo_,g);
}

template<class V>
void BoundConstraint<V>::pruneActive( V& v, const V& x, Real eps ) const {
  if (isActivated()) {
    auto epsn{std::min(scale_*eps,min_diff_)}; 
    auto f = [epsn]( auto ve, auto xe, auto ue, auto le ) { 
      return ( ((ue-xe)<=epsn) || ((xe-le)<=epsn) ) ? 0 : ve;
    };
    Elementwise::eval_function(v,f,x,*x_up_,*x_lo_);
  }
}

template<class V>
void BoundConstraint<V>::pruneActive( V& v, const dual_t<V>& g, const V& x, Real eps ) const {
  if (isActivated()) {
    auto epsn{std::min(scale_*eps,min_diff_)}; 
    auto f = [epsn]( auto ve, auto xe, auto ue, auto le, auto ge ) { 
      return ( ( ((ue-xe)<=epsn) && (ge<0) ) || ( ((xe-le)<=epsn) && (ge>0) ) ) ? 0 : ve;
    };
    Elementwise::eval_function(v,f,x,*x_up_,*x_lo_,g);
  }
}


template<class V>
const V& BoundConstraint<V>::getLowerBound( void ) const {
  return *x_lo_;
}

template<class V>
const V& BoundConstraint<V>::getUpperBound( void ) const {
  return *x_up_;
}

template<class V>
bool BoundConstraint<V>::isFeasible( const V& v ) const {
  auto is_between = []( auto l, auto v, auto u ) {
    return ( (l<v) && (v<u) ) ? Real(1) : Real(0);
  }; 
  auto all_true = make_product(v);
  return (Elementwise::eval_function_and_reduce(all_true,is_between,*x_lo_,v,*x_up_) != 0);
}

template<class V>
void BoundConstraint<V>::activateLower(void) {
  Lactivated_ = true;
}

template<class V>
void BoundConstraint<V>::activateUpper(void) {
  Uactivated_ = true;
}

template<class V>
void BoundConstraint<V>::activate(void) {
  activateLower();
  activateUpper();
}

template<class V>
void BoundConstraint<V>::deactivateLower(void) {
  Lactivated_ = false;
}

template<class V>
void BoundConstraint<V>::deactivateUpper(void) {
  Uactivated_ = false;
}

template<class V>
void BoundConstraint<V>::deactivate(void) {
  deactivateLower();
  deactivateUpper();
}

template<class V>
bool BoundConstraint<V>::isLowerActivated(void) const {
  return Lactivated_;
}

template<class V>
bool BoundConstraint<V>::isUpperActivated(void) const {
  return Uactivated_;
}

template<class V>
bool BoundConstraint<V>::isActivated(void) const {
  return (isLowerActivated() || isUpperActivated());
}


template<class V>
void BoundConstraint<V>::pruneLowerInactive( V& v, const V& x, Real eps ) {
  if( !isLowerActivated() ) {
    auto epsn{std::min(scale_*eps,min_diff_)}; 
    auto f = [epsn]( auto ve, auto xe, auto le) { 
      return ( (xe-le)<=epsn ) ? ve : 0; 
    };
    Elementwise::eval_function(v,f,x,*x_lo_);
  }
}

template<class V>
void BoundConstraint<V>::pruneLowerInactive( V& v, const dual_t<V>& g, const V& x, Real eps ) {
  if( !isLowerActivated() ) {
    auto epsn{std::min(scale_*eps,min_diff_)}; 
    auto f = [epsn]( auto ve, auto xe, auto le, auto ge ) { 
      return ( ((xe-le)<=epsn) && ge>0 ) ? ve : 0; 
    };
    Elementwise::eval_function(v,f,x,*x_lo_,g);
  }
}

template<class V>
void BoundConstraint<V>::pruneUpperInactive( V& v, const V& x, Real eps ) {
  auto epsn{std::min(scale_*eps,min_diff_)}; 
  auto f = [epsn]( auto ve, auto xe, auto ue) { 
    return ( (ue-xe)<=epsn ) ? ve : 0; 
  };
  Elementwise::eval_function(v,f,x,*x_up_);
}

template<class V>
void BoundConstraint<V>::pruneUpperInactive( V& v, const dual_t<V>& g, const V& x, Real eps ) {
  auto epsn{std::min(scale_*eps,min_diff_)}; 
  auto f = [epsn]( auto ve, auto xe, auto ue, auto ge ) { 
    return ( ((ue-xe)<=epsn) && ge<0 ) ? ve : 0; 
  };
  Elementwise::eval_function(v,f,x,*x_up_,g);
}

template<class V>
void BoundConstraint<V>::pruneInactive( V& v, const V& x, Real eps ) const {
  if (isActivated()) {
    auto epsn{std::min(scale_*eps,min_diff_)}; 
    auto f = [epsn]( auto ve, auto xe, auto ue, auto le ) { 
      return ( ((ue-xe)<=epsn) || ((xe-le)<=epsn) ) ? ve : 0;
    };
    Elementwise::eval_function(v,f,x,*x_up_,*x_lo_);
  }
}

template<class V>
void BoundConstraint<V>::pruneInactive( V& v, const dual_t<V>& g, const V& x, Real eps ) const {
  if (isActivated()) {
    auto epsn{std::min(scale_*eps,min_diff_)}; 
    auto f = [epsn]( auto ve, auto xe, auto ue, auto le, auto ge ) { 
      return ( ( ((ue-xe)<=epsn) && (ge<0) ) || ( ((xe-le)<=epsn) && (ge>0) ) ) ? ve : 0;
    };
    Elementwise::eval_function(v,f,x,*x_up_,*x_lo_,g);
  }
}


} // namespace XROL

