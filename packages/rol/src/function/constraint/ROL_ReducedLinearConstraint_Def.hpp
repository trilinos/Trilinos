// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_REDUCED_LINEAR_CONSTRAINT_DEF_H
#define ROL_REDUCED_LINEAR_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
ReducedLinearConstraint<Real>::ReducedLinearConstraint(const Ptr<Constraint<Real>> &con,
                                                       const Ptr<BoundConstraint<Real>> &bnd,
                                                       const Ptr<const Vector<Real>> &x)
  : con_(con), bnd_(bnd), x_(x), prim_(x->clone()) {}

template<typename Real>
void ReducedLinearConstraint<Real>::setX(const Ptr<const Vector<Real>> &x) {
  x_ = x;
}

template<typename Real>
void ReducedLinearConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
  const Real zero(0);
  prim_->set(x);
  bnd_->pruneActive(*prim_,*x_,zero);
  con_->value(c,*prim_,tol);
}

template<typename Real>
void ReducedLinearConstraint<Real>::applyJacobian(Vector<Real> &jv,
                                            const Vector<Real> &v,
                                            const Vector<Real> &x, Real &tol) {
  const Real zero(0);
  prim_->set(v);
  bnd_->pruneActive(*prim_,*x_,zero);
  con_->applyJacobian(jv,*prim_,x,tol);
}

template<typename Real>
void ReducedLinearConstraint<Real>::applyAdjointJacobian(Vector<Real> &jv,
                                                   const Vector<Real> &v,
                                                   const Vector<Real> &x, Real &tol) {
  const Real zero(0);
  con_->applyAdjointJacobian(jv,v,x,tol);
  bnd_->pruneActive(jv,*x_,zero);
}

template<typename Real>
void ReducedLinearConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv,
                                                  const Vector<Real> &u,
                                                  const Vector<Real> &v,
                                                  const Vector<Real> &x, Real &tol) {
  ahuv.zero();
}

} // namespace ROL

#endif // ROL_REDUCED_LINEAR_CONSTRAINT_H
