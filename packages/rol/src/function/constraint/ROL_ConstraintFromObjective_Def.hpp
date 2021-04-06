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

#ifndef ROL_CONSTRAINTFROMOBJECTIVE_DEF_H
#define ROL_CONSTRAINTFROMOBJECTIVE_DEF_H

namespace ROL {

template<typename Real> 
ConstraintFromObjective<Real>::ConstraintFromObjective( const Ptr<Objective<Real>> &obj, const Real offset ) :
  obj_(obj), dualVector_(nullPtr), offset_(offset), isDualInitialized_(false) {}

template<typename Real> 
const Ptr<Objective<Real>> ConstraintFromObjective<Real>::getObjective(void) const { return obj_; }

template<typename Real> 
void ConstraintFromObjective<Real>::setParameter( const std::vector<Real> &param ) {
  obj_->setParameter(param);
  Constraint<Real>::setParameter(param);
}

template<typename Real> 
void ConstraintFromObjective<Real>::update( const Vector<Real>& x, UpdateType type, int iter ) {
  obj_->update(x,type,iter);
}

template<typename Real> 
void ConstraintFromObjective<Real>::update( const Vector<Real>& x, bool flag, int iter ) {
  obj_->update(x,flag,iter);
}

template<typename Real> 
void ConstraintFromObjective<Real>::value( Vector<Real>& c, const Vector<Real>& x, Real& tol ) {
  setValue(c, obj_->value(x,tol) - offset_ ); 
}

template<typename Real> 
void ConstraintFromObjective<Real>::applyJacobian( Vector<Real>& jv, const Vector<Real>& v, const Vector<Real>& x, Real& tol ) {
  if ( !isDualInitialized_ ) {
    dualVector_ = x.dual().clone();
    isDualInitialized_ = true;
  }
  obj_->gradient(*dualVector_,x,tol);
  //setValue(jv,v.dot(dualVector_->dual()));
  setValue(jv,v.apply(*dualVector_));
}

template<typename Real> 
void ConstraintFromObjective<Real>::applyAdjointJacobian( Vector<Real>& ajv, const Vector<Real>& v, const Vector<Real>& x, Real& tol ) {
  obj_->gradient(ajv,x,tol);
  ajv.scale(getValue(v));
}

template<typename Real> 
void ConstraintFromObjective<Real>::applyAdjointHessian( Vector<Real>& ahuv, const Vector<Real>& u, const Vector<Real>& v, const Vector<Real>& x, Real& tol ) {
  obj_->hessVec(ahuv,v,x,tol);
  ahuv.scale(getValue(u));
}

template<typename Real> 
Real ConstraintFromObjective<Real>::getValue( const Vector<Real>& x ) { 
  return dynamic_cast<const SingletonVector<Real>&>(x).getValue(); 
}

template<typename Real> 
void ConstraintFromObjective<Real>::setValue( Vector<Real>& x, Real val ) {
  dynamic_cast<SingletonVector<Real>&>(x).setValue(val);
}

} // namespace ROL

#endif // ROL_CONSTRAINTFROMOBJECTIVE_H
