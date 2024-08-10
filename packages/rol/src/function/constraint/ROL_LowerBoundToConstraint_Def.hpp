// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LOWER_BOUND_TO_CONSTRAINT_DEF_H
#define ROL_LOWER_BOUND_TO_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
LowerBoundToConstraint<Real>::LowerBoundToConstraint(BoundConstraint<Real> &bnd) {
  lo_ = bnd.getLowerBound()->clone();
  lo_->set(*bnd.getLowerBound());
}

template<typename Real>
LowerBoundToConstraint<Real>::LowerBoundToConstraint(const Vector<Real> &lo) {
  lo_ = lo.clone();
  lo_->set(lo);
}

template<typename Real>
void LowerBoundToConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
  const Real one(1);
  c.set(x);
  c.axpy(-one,*lo_);
}

template<typename Real>
void LowerBoundToConstraint<Real>::applyJacobian(Vector<Real> &jv,
                                                 const Vector<Real> &v,
                                                 const Vector<Real> &x,
                                                 Real &tol) {
  jv.set(v);
}

template<typename Real>
void LowerBoundToConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,
                                                        const Vector<Real> &v,
                                                        const Vector<Real> &x,
                                                        Real &tol) {
  ajv.set(v);
}

template<typename Real>
void LowerBoundToConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv,
                                                       const Vector<Real> &u,
                                                       const Vector<Real> &v,
                                                       const Vector<Real> &x,
                                                       Real &tol) {
  ahuv.zero();
}

}

#endif
