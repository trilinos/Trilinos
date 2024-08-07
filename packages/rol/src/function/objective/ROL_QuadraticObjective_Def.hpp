// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_QUADRATIC_OBJECTIVE_DEF_H
#define ROL_QUADRATIC_OBJECTIVE_DEF_H

namespace ROL {

template<typename Real>
QuadraticObjective<Real>::QuadraticObjective(const Ptr<const LinearOperator<Real>> &H,
                                             const Ptr<const Vector<Real>>         &g,
                                             Real                                   c)
  : H_(H), g_(g), c_(c) {
  tmp_ = g_->clone();
}

template<typename Real>
Real QuadraticObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  H_->apply(*tmp_,x,tol);
  tmp_->scale(static_cast<Real>(0.5));
  tmp_->plus(*g_);
  //return x.dot(tmp_->dual()) + c_;
  return x.apply(*tmp_) + c_;
}

template<typename Real>
void QuadraticObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  H_->apply(g,x,tol);
  g.plus(*g_);
}

template<typename Real>
void QuadraticObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  H_->apply(hv,v,tol);
}

template<typename Real>
void QuadraticObjective<Real>::invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  H_->applyInverse(hv,v,tol);
}

} // namespace ROL

#endif
