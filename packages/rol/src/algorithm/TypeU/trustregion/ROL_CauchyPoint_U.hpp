// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CAUCHYPOINT_U_H
#define ROL_CAUCHYPOINT_U_H

/** \class ROL::CauchyPoint_U
    \brief Provides interface for the Cauchy point trust-region subproblem solver.
*/

#include "ROL_TrustRegion_U.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class CauchyPoint_U : public TrustRegion_U<Real> {
private:

  Ptr<Vector<Real>> dual_;

public:

  CauchyPoint_U() {}

  void initialize(const Vector<Real> &x, const Vector<Real> &g) {
    dual_ = g.clone();
  }

  void solve( Vector<Real> &s, Real &snorm, Real &pRed,
              int &iflag, int &iter, const Real del,
              TrustRegionModel_U<Real> &model) {
    const Real zero(0), half(0.5);
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    // Set step to (projected) gradient
    s.set(model.getGradient()->dual());
    // Apply (reduced) Hessian to (projected) gradient
    model.hessVec(*dual_,s,s,tol);
    Real gnorm  = s.norm();
    Real gnorm2 = gnorm*gnorm;
    //Real gBg    = dual_->dot(s.dual());
    Real gBg    = dual_->apply(s);
    Real alpha  = gnorm2/gBg;
    if ( alpha*gnorm >= del || gBg <= zero ) {
      alpha = del/gnorm;
    }
    s.scale(-alpha);
    snorm = alpha*gnorm;
    iflag = 0;
    iter  = 0;
    pRed  = alpha*(gnorm2 - half*alpha*gBg);
  }
};

} // namespace ROL

#endif
