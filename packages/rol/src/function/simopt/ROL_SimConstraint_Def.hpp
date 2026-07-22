// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINT_STATE_DEF_H
#define ROL_CONSTRAINT_STATE_DEF_H

namespace ROL {

template<typename Real>
SimConstraint<Real>::SimConstraint(const Ptr<Constraint_SimOpt<Real>> &con,
                                   const Ptr<const Vector<Real>> &z,
                                   bool inSolve) : con_(con), z_(z), inSolve_(inSolve), init_(false) {}

template<typename Real>
void SimConstraint<Real>::update( const Vector<Real> &u, bool flag, int iter ) {
  con_->update_1(u,flag,iter);
  //con_->update(u,*z_,flag,iter);
}

template<typename Real>
void SimConstraint<Real>::update( const Vector<Real> &u, UpdateType type, int iter ) {
  if (inSolve_) con_->solve_update(u,*z_,type,iter);
  else          con_->update_1(u,type,iter);
}

template<typename Real>
void SimConstraint<Real>::value(Vector<Real> &c,const Vector<Real> &u,Real &tol) {
  con_->value(c,u,*z_,tol);
}

template<typename Real>
void SimConstraint<Real>::applyJacobian(Vector<Real> &jv,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
  con_->applyJacobian_1(jv,v,u,*z_,tol);
}

template<typename Real>
void SimConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
  con_->applyAdjointJacobian_1(ajv,v,u,*z_,tol);
}

template<typename Real>
void SimConstraint<Real>::applyAdjointHessian(Vector<Real> &ahwv,const Vector<Real> &w,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
  con_->applyAdjointHessian_11(ahwv,w,v,u,*z_,tol);
}

template<typename Real>
void SimConstraint<Real>::applyPreconditioner(Vector<Real> &pv,const Vector<Real> &v,const Vector<Real> &u,const Vector<Real> &g,Real &tol) {
  if (!init_) {
    ijv_ = u.clone();
    init_ = true;
  }
  con_->applyInverseJacobian_1(*ijv_,v,u,*z_,tol);
  con_->applyInverseAdjointJacobian_1(pv,ijv_->dual(),u,*z_,tol);
}

template<typename Real>
void SimConstraint<Real>::setParameter(const std::vector<Real> &param) {
  con_->setParameter(param);
  Constraint<Real>::setParameter(param);
}

} // namespace ROL

#endif
