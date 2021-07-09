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

#ifndef ROL_OBJECTIVE_FROM_CONSTRAINT_DEF_H
#define ROL_OBJECTIVE_FROM_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
ObjectiveFromConstraint<Real>::ObjectiveFromConstraint( const Ptr<Constraint<Real>> &con, 
                                                        const Vector<Real> &l ) :
  con_(con), l_(l.clone()), c_(l.dual().clone()) {
  l_->set(l);
}

template<typename Real>
void ObjectiveFromConstraint<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  con_->update(x,type,iter);
}

template<typename Real>
void ObjectiveFromConstraint<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  con_->update(x,flag,iter);
}

template<typename Real>
Real ObjectiveFromConstraint<Real>::value( const Vector<Real> &x, Real &tol ) {
  con_->value(*c_,x,tol);
  //return l_->dot(c_->dual());  
  return l_->apply(*c_);  
}

template<typename Real>
void ObjectiveFromConstraint<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  con_->applyAdjointJacobian(g,*l_,x,tol);
}

template<typename Real>
void ObjectiveFromConstraint<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  con_->applyAdjointHessian(hv,*l_,v,x,tol);
}

template<typename Real>
void ObjectiveFromConstraint<Real>::updateMultiplier( const Vector<Real> &l ) {
  l_->set(l);
}

} // namespace ROL

#endif // ROL_OBJECTIVE_FROM_CONSTRAINT_DEF_H
