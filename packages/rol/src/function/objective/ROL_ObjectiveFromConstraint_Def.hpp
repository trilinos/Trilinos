// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
