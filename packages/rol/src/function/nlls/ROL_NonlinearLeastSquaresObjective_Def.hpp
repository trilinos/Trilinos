// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_NONLINEARLEASTSQUARESOBJECTIVE_DEF_H
#define ROL_NONLINEARLEASTSQUARESOBJECTIVE_DEF_H

namespace ROL {

template<typename Real>
NonlinearLeastSquaresObjective<Real>::NonlinearLeastSquaresObjective(const Ptr<Constraint<Real>> &con,
                                                                     const Vector<Real> &optvec,
                                                                     const Vector<Real> &convec,
                                                                     const bool GNH)
  : con_(con), GaussNewtonHessian_(GNH) {
  c1_ = convec.clone(); c1dual_ = c1_->dual().clone();
  c2_ = convec.clone();
  x_  = optvec.dual().clone();
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  con_->update(x,type,iter);
  con_->value(*c1_,x,tol);
  c1dual_->set(c1_->dual());
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  con_->update(x,flag,iter);
  con_->value(*c1_,x,tol);
  c1dual_->set(c1_->dual());
}

template<typename Real>
Real NonlinearLeastSquaresObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  Real half(0.5);
  return half*(c1_->dot(*c1_));
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  con_->applyAdjointJacobian(g,*c1dual_,x,tol);
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  con_->applyJacobian(*c2_,v,x,tol);
  con_->applyAdjointJacobian(hv,c2_->dual(),x,tol);
  if ( !GaussNewtonHessian_ ) {
    con_->applyAdjointHessian(*x_,*c1dual_,v,x,tol);
    hv.plus(*x_);
  }
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  con_->applyPreconditioner(Pv,v,x,x.dual(),tol);
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::setParameter(const std::vector<Real> &param) {
  Objective<Real>::setParameter(param);
  con_->setParameter(param);
}

} // namespace ROL

#endif
