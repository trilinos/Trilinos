
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
#include "XROL_Vector.hpp"


namespace XROL {
template<class X>
magnitude_t<X>
Objective<X>::dirDeriv( const X& x, const X& d, magnitude_t<X>& tol ) {
  vec_->allocate_x(x);
  X& xd = vec_->x();
  set(xd,x);
  axpy(xd,tol,d);
  return ( value(xd,tol) - value(x,tol) )/tol;
} // dirDeriv


template<class X>
void Objective<X>::gradient( dual_t<X>& g, const X& x, magnitude_t<X>& tol ) {
  vec_->allocate_x(x);
  vec_->allocate_g(x);
  X& exi = vec_->x();
  X& egi = vec_->g();
  fill(g,0);
  magnitude_t<X> deriv = 0.0;
  magnitude_t<X> h     = 0.0;   
  magnitude_t<X> xi    = 0.0;
  auto eps = std::numeric_limits<magnitude_t<X>>::epsilon();
  for( index_t<X> i=0; i<dimension(g); ++i ) {
    basis(exi,i);
    basis(egi,i);
    xi     = std::abs(dot(x,exi));
    h      = ( xi < eps ? 1.0 : xi )*tol;
    deriv  = dirDeriv(x,exi,h);
    axpy(g,deriv,egi);
  }
} // gradient


template<class X>
void Objective<X>::hessVec( dual_t<X> &hv, const X &v, const X &x, magnitude_t<X>& tol ) {

  vec_->allocate_g(x);
  vec_->allocate_g(hv);
  magnitude_t<X> _one(1.0);
  magnitude_t<X> _zero(0.0);
  auto eps = std::numeric_limits<magnitude_t<X>>::epsilon();
   
  if( norm(v) == _zero ) {
    fill(hv,_zero);
  }
  else {
    auto gtol = std::sqrt(eps);
    auto h = std::max(_one,norm(x)/norm(v))*tol;
    dual_t<X>& g = vec_->g();
    gradient(g,x,gtol);


    X& xnew = vec_->x();
    set(xnew,x);
    axpy(xnew,h,v);
    update(xnew);

    fill(hv,_zero);
    gradient(hv,xnew,gtol);
    axpy(hv,-_one,g);
    scale(hv,_one/h);
  }
} // hessVec


template<class X>
void Objective<X>::update( const X& x ) {
  ignore(x);
}  // update


template<class X>
void Objective<X>::invHessVec( X& hv, const dual_t<X>& v, const X& x, magnitude_t<X>& tol ) {
  ignore(hv,v,x,tol);
} // invHessVec

template<class X>
void Objective<X>::precond( X& Pv, const dual_t<X>& v, const X& x, magnitude_t<X>& tol ) {
  ignore(Pv,v,x,tol);
} // precond


} // namespace XROL
