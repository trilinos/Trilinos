// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALED_OBJECTIVE_DEF_HPP
#define ROL_SCALED_OBJECTIVE_DEF_HPP

namespace ROL {

template<typename Real>
void ScaledObjective<Real>::update(const Vector<Real> &x, UpdateType type, int iter) {
  obj_->update(x,type,iter);
}

template<typename Real>
void ScaledObjective<Real>::setParameter(const std::vector<Real> &param) {
  Objective<Real>::setParameter(param);
  obj_->setParameter(param);
}

template<typename Real>
Real ScaledObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  return scale_ * obj_->value(x,tol);
}

template<typename Real>
void ScaledObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  obj_->gradient(g,x,tol);
  g.scale(scale_);
}

template<typename Real>
void ScaledObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  obj_->hessVec(hv,v,x,tol);
  hv.scale(scale_);
}

template<typename Real>
void ScaledObjective<Real>::invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  obj_->invHessVec(hv,v,x,tol);
  hv.scale(static_cast<Real>(1)/scale_);
}

template<typename Real>
void ScaledObjective<Real>::precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  obj_->precond(Pv,v,x,tol);
  Pv.scale(static_cast<Real>(1)/scale_);
}

} // End ROL Namespace

#endif
