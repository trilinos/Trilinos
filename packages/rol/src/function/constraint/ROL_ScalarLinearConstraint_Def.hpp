// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AFFINE_HYPERPLANE_EQUALITY_CONSTRAINT_DEF_H
#define ROL_AFFINE_HYPERPLANE_EQUALITY_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
ScalarLinearConstraint<Real>::ScalarLinearConstraint(const Ptr<const Vector<Real>> &a,
                                                     const Real b)
  : a_(a), b_(b) {}

template<typename Real>
void ScalarLinearConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
  SingletonVector<Real> &cc = dynamic_cast<SingletonVector<Real>&>(c);
  //cc.setValue(a_->dot(x.dual()) - b_);
  cc.setValue(a_->apply(x) - b_);
}

template<typename Real>
void ScalarLinearConstraint<Real>::applyJacobian(Vector<Real> &jv,
                                           const Vector<Real> &v,
                                           const Vector<Real> &x, Real &tol) {
  SingletonVector<Real> &jc = dynamic_cast<SingletonVector<Real>&>(jv);
  //jc.setValue(a_->dot(v.dual()));
  jc.setValue(a_->apply(v));
}

template<typename Real>
void ScalarLinearConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,
                                                  const Vector<Real> &v,
                                                  const Vector<Real> &x, Real &tol) {
  const SingletonVector<Real>&    vc = dynamic_cast<const SingletonVector<Real>&>(v);
  ajv.set(*a_);
  ajv.scale(vc.getValue());
}

template<typename Real>
void ScalarLinearConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv,
                                                 const Vector<Real> &u,
                                                 const Vector<Real> &v,
                                                 const Vector<Real> &x, Real &tol) {
  ahuv.zero();
}

template<typename Real>
std::vector<Real> ScalarLinearConstraint<Real>::solveAugmentedSystem(Vector<Real> &v1,
                                                                     Vector<Real> &v2,
                                                               const Vector<Real> &b1,
                                                               const Vector<Real> &b2,
                                                               const Vector<Real> &x, Real &tol) {
  SingletonVector<Real>&       v2c = dynamic_cast<SingletonVector<Real>&>(v2);
  const SingletonVector<Real>& b2c = dynamic_cast<const SingletonVector<Real>&>(b2);

  //v2c.setValue( (a_->dot(b1.dual()) - b2c.getValue() )/a_->dot(*a_) );
  v2c.setValue( (a_->apply(b1) - b2c.getValue() )/a_->dot(*a_) );
  v1.set(b1.dual());
  v1.axpy(-v2c.getValue(),a_->dual());

  std::vector<Real> out;
  return out;
}

} // namespace ROL

#endif
