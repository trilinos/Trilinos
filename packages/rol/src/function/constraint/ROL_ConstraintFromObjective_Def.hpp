// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
