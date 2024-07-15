// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OBJECTIVE_FSSOLVER_DEF_H
#define ROL_OBJECTIVE_FSSOLVER_DEF_H

namespace ROL {

template<typename Real>
Real Objective_FSsolver<Real>::value( const Vector<Real> &u, Real &tol ) {
  return static_cast<Real>(0.5)*u.dot(u);
}

template<typename Real>
void Objective_FSsolver<Real>::gradient( Vector<Real> &g, const Vector<Real> &u, Real &tol ) {
  g.set(u.dual());
}

template<typename Real>
void Objective_FSsolver<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, Real &tol ) {
  hv.set(v.dual());
}

} // namespace ROL

#endif
