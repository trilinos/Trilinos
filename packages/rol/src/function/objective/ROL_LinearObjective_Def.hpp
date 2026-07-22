// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINEAR_OBJECTIVE_DEF_H
#define ROL_LINEAR_OBJECTIVE_DEF_H

namespace ROL {

template<typename Real>
LinearObjective<Real>::LinearObjective(const Ptr<const Vector<Real>> &cost) : cost_(cost) {
  dual_cost_ = cost_->dual().clone();
  dual_cost_->set(cost_->dual());
}

template<typename Real>
Real LinearObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  return x.dot(*dual_cost_);
}

template<typename Real>
void LinearObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  g.set(*cost_);
}

template<typename Real>
void LinearObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  hv.zero();
}

} // namespace ROL

#endif
