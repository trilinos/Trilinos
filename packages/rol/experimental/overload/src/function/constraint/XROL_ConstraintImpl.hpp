
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

#include "XROL.hpp"


namespace XROL {

template<class X, class C>
void Constraint<X,C>::Constraint() : 
  activated_(true), c_allocated_(false) {}

template<class X, class C>
void Constraint<X,C>::~Constraint(){}

template<class X, class C>
void Constraint<X,C>::update( const X& x ) { ignore(x); }

template<class X, class C>
void Constraint<X,C>::applyJacobian( C& jv,
                                     const X& v,
                                     const X& x, 
                                     magnitude_t<X>& tol ) {
  magnitude_t<X> one(1);
  auto eps = std::numeric_limits<magnitude_t<X>>::epsilon();
  magnitude_t<X> ctol = std::sqrt(eps);

  auto h = std::max(one,norm(x)/norm(v))*tol;

  convec_->allocate_x(x,2);
  convec_->allocate_c(jv,3);

  C& c = convec_->c(0);
  value(c,x,ctol);

  X& xnew = concec_->x(0);
  set(xnew,x);
  axpy(xnew,h,v);
  update(xnew);

  fill(jv,0);
  value(jv,xnew,ctol);

  axpy(jv,-one,c);
  scale(jv,one/h);
}



template<class X, class C> 
void Constraint<X,C>::applyAdjointJacobian( dual_t<X>& ajv,  
                                            const dual_t<C>& v,
                                            const X& x, 
                                            const C& c,
                                            magnitude_t<X>& tol ) {
  
  magnitude_t<X> one(1);
  magnitude_t<X> h(0);
  auto eps = std::numeric_limits<magnitude_t<X>>::epsilon();
  magnitude_t<X> ctol = std::sqrt(eps);

  convec_->allocate_x(x,2);
  convec_->allocate_g(ajv);
  convec_->allocate_c(jv,3);
 
  X& xnew  = convec_->x(0);
  X& ex    = convec_->x(1);  
  dual_t<X>& eajv = convec_->g();
  C& cnew  = convec_->c(0);
  C& c0    = convec_->c(1);
  C& vdual = convec_->c(2);

  dual(vdual,v);

  magnitude_t<X> one(1);
  magnitude_t<X> h(0);
  magnitude_t<X> ctol; 

  value(c0,x,ctol);

  fill(ajv,0);

  for( index_t<X> i=0; i<ajv.dimension(); ++i ) {
    basis(ex,i);
    basis(eajv,i);
    h = std::max(one,norm(x)/norm(ex))*tol;
    set(xnew,x);
    update(xnew);
    value(cnew,xnew,ctol);
    axpy(cnew,-one,c0);
    scale(one/h);
    axpy(ajv,dot(cnew,vdual),eajv);
  }

} 

template<class X, class C>
void Constraint<X,C>::applyAdjointHessian( dual_t<X>& ahuv,
                                           const dual_t<C>& u,
                                           const X& v,
                                           const X& x, 
                                           magnitude_t<X>& tol ) {
  magnitude_t<X> one(1);
  magnitude_t<X> h = std::max(one,norm(x)/norm(v))*tol;
 
  convec_->allocate_x(x,2);
  convec_->allocate_g(ajv);

  dual_t<X>& aju = convec_->g();
  applyAdjointJacobian(aju,u,x,tol);
  X& xnew = convec_->x();
  set(xnew,x);
  axpy(xnew,h,v);
  update(xnew);

}

/*
std::vector<magnitude_t<X>> 
Constraint<X,C>::solveAugmentedSystem( X&                v1,
                                       dual_t<C>&        v2,
                                       const dual_t<C>&  b1,
                                       const C&          b2,
                                       const X&          x,
                                       magnitude_t<X>&   tol ) {

  magnitude_t<X> zero(0), one(1);
  
}
*/


template<class X, class C>
void Constraint<X,C>::applyPreconditioner( dual_t<C>& pv, 
                                           const C& v, 
                                           const X& x, 
                                           const dual_t<X>& g,
                                           magnitude_t<X>& tol ) {
  ignore(pv,v,x,g,tol);
}


template<class X, class C>
void Constraint<X,C>::activate( void ) {
  activated_ = true;
}

template<class X, class C>
void Constraint<X,C>::deactivate( void ) {
  activated_ = false;
}

} // namespace XROL

