// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINEARCONSTRAINT_DEF_H
#define ROL_LINEARCONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
LinearConstraint<Real>::LinearConstraint(const Ptr<const LinearOperator<Real>> &A,
                                         const Ptr<const Vector<Real>> &b) : A_(A), b_(b) {}

template<typename Real>
void LinearConstraint<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {}

template<typename Real>
void LinearConstraint<Real>::update( const Vector<Real> &x, bool flag, int iter ) {}

template<typename Real>
void LinearConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
  A_->apply(c,x,tol);
  c.plus(*b_);
}

template<typename Real>
void LinearConstraint<Real>::applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  A_->apply(jv,v,tol);
}

template<typename Real>
void LinearConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  A_->applyAdjoint(ajv,v,tol);
}

template<typename Real>
void LinearConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &dualv, Real &tol) {
  A_->applyAdjoint(ajv,v,tol);
}

template<typename Real>
void LinearConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  ahuv.zero();
}

template<typename Real>
Ptr<Vector<Real>> LinearConstraint<Real>::createRangeSpaceVector(void) const {
  return b_->clone();
}

} // namespace ROL

#endif
