// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PQNOBJECTIVEDEF_H
#define ROL_PQNOBJECTIVEDEF_H


namespace ROL {

template<typename Real>
PQNObjective<Real>::PQNObjective(const Ptr<Secant<Real>> &secant,
                                 const Vector<Real> &x,
                                 const Vector<Real> &g)
  : secant_(secant), x_(x.clone()), g_(g.clone()), pwa_(x.clone()), dwa_(g.clone()) {
  setAnchor(x,g);
}

template<typename Real>
Real PQNObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  pwa_->set(x);
  pwa_->axpy(static_cast<Real>(-1),*x_);
  secant_->applyB(*dwa_, *pwa_);
  dwa_->scale(static_cast<Real>(0.5));
  dwa_->plus(*g_);
  return dwa_->apply(*pwa_);
}

template<typename Real>
void PQNObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  pwa_->set(x);
  pwa_->axpy(static_cast<Real>(-1),*x_);
  secant_->applyB(g, *pwa_);
  g.plus(*g_);
}

template<typename Real>
void PQNObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  secant_->applyB(hv, v);
}

template<typename Real>
void PQNObjective<Real>::setAnchor(const Vector<Real> &x, const Vector<Real> &g) {
  x_->set(x); g_->set(g);
}

} // namespace ROL

#endif
