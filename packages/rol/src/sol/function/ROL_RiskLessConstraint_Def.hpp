// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RISKLESS_CONSTRAINT_DEF_H
#define ROL_RISKLESS_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
RiskLessConstraint<Real>::RiskLessConstraint(const Ptr<Constraint<Real>> &con)
  : con_(con) {}

template<typename Real>
void RiskLessConstraint<Real>::update(const Vector<Real> &x, UpdateType type, int iter) {
  Ptr<const Vector<Real>> x0 = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  con_->update(*x0,type,iter);
}

template<typename Real>
void RiskLessConstraint<Real>::update(const Vector<Real> &x, bool flag, int iter) {
  Ptr<const Vector<Real>> x0 = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  con_->update(*x0,flag,iter);
}

template<typename Real>
void RiskLessConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
  Ptr<const Vector<Real>> x0 = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  con_->value(c,*x0,tol);
}

template<typename Real>
void RiskLessConstraint<Real>::applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  Ptr<const Vector<Real>> x0 = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  Ptr<const Vector<Real>> v0 = dynamic_cast<const RiskVector<Real>&>(v).getVector();
  con_->applyJacobian(jv,*v0,*x0,tol);
}

template<typename Real>
void RiskLessConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  Ptr<const Vector<Real>> x0 = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  Ptr<Vector<Real>> ajv0 = dynamic_cast<RiskVector<Real>&>(ajv).getVector();
  con_->applyAdjointJacobian(*ajv0,v,*x0,tol);
}

template<typename Real>
void RiskLessConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  Ptr<const Vector<Real>> x0 = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  Ptr<const Vector<Real>> v0 = dynamic_cast<const RiskVector<Real>&>(v).getVector();
  Ptr<Vector<Real>> ahuv0 = dynamic_cast<RiskVector<Real>&>(ahuv).getVector();
  con_->applyAdjointHessian(*ahuv0,u,*v0,*x0,tol);
}

template<typename Real>
void RiskLessConstraint<Real>::setParameter(const std::vector<Real> &param) {
  Constraint<Real>::setParameter(param);
  con_->setParameter(param);
}

} // namespace ROL

#endif
