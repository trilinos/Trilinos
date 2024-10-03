// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BOUND_TO_CONSTRAINT_H
#define ROL_BOUND_TO_CONSTRAINT_H

namespace ROL {

template<typename Real>
BoundToConstraint<Real>::BoundToConstraint(BoundConstraint<Real> &bnd) {
  lo_ = makePtr<LowerBoundToConstraint<Real>>(bnd);
  up_ = makePtr<UpperBoundToConstraint<Real>>(bnd);
  tmp_ = x.clone();
}

template<typename Real>
BoundToConstraint<Real>::BoundToConstraint(const Vector<Real> &lo, const Vector<Real> &up) {
  lo_ = makePtr<LowerBoundToConstraint<Real>>(lo);
  up_ = makePtr<UpperBoundToConstraint<Real>>(up);
  tmp_ = lo.clone();
}

template<typename Real>
void BoundToConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
  Vector<Real> &c0 = *(dynamic_cast<PartitionedVector<Real>&>(c).get(0));
  Vector<Real> &c1 = *(dynamic_cast<PartitionedVector<Real>&>(c).get(1));
  lo_->value(c0,x,tol);
  up_->value(c1,x,tol);
}

template<typename Real>
void BoundToConstraint<Real>::applyJacobian(Vector<Real> &jv,
                   const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  Vector<Real> &jv0 = *(dynamic_cast<PartitionedVector<Real>&>(jv).get(0));
  Vector<Real> &jv1 = *(dynamic_cast<PartitionedVector<Real>&>(jv).get(1));
  lo_->applyJacobian(jv0,v,x,tol);
  up_->applyJacobian(jv1,v,x,tol);
}

template<typename Real>
void BoundToConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,
                   const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  const Vector<Real> &v0 = *(dynamic_cast<const PartitionedVector<Real>&>(v).get(0));
  const Vector<Real> &v1 = *(dynamic_cast<const PartitionedVector<Real>&>(v).get(1));
  lo_->applyAdjointJacobian(ajv,v0,x,tol);
  up_->applyAdjointJacobian(*tmp_,v1,x,tol);
  ajv.plus(*tmp_); 
}

template<typename Real>
void BoundToConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv,
                   const Vector<Real> &u, const Vector<Real> &v,
                   const Vector<Real> &x, Real &tol) {
  ahuv.zero();
}

}

#endif
