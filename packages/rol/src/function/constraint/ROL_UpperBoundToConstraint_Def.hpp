// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_UPPER_BOUND_TO_CONSTRAINT_DEF_H
#define ROL_UPPER_BOUND_TO_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
UpperBoundToConstraint<Real>::UpperBoundToConstraint(BoundConstraint<Real> &bnd) {
  up_ = bnd.getUpperBound()->clone();
  up_->set(*bnd.getUpperBound());
}

template<typename Real>
UpperBoundToConstraint<Real>::UpperBoundToConstraint(const Vector<Real> &up) {
  up_ = up.clone();
  up_->set(up);
}

template<typename Real>
void UpperBoundToConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
  const Real one(1);
  c.set(*up_);
  c.axpy(-one,x);
}

template<typename Real>
void UpperBoundToConstraint<Real>::applyJacobian(Vector<Real> &jv,
                                                 const Vector<Real> &v,
                                                 const Vector<Real> &x,
                                                 Real &tol) {
  const Real one(1);
  jv.set(v);
  jv.scale(-one);
}

template<typename Real>
void UpperBoundToConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,
                                                        const Vector<Real> &v,
                                                        const Vector<Real> &x,
                                                        Real &tol) {
  const Real one(1);
  ajv.set(v);
  ajv.scale(-one);
}

template<typename Real>
void UpperBoundToConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv,
                                                       const Vector<Real> &u,
                                                       const Vector<Real> &v,
                                                       const Vector<Real> &x,
                                                       Real &tol) {
  ahuv.zero();
}

}

#endif
