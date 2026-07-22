// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_NEWTON_U_H
#define ROL_NEWTON_U_H

#include "ROL_DescentDirection_U.hpp"
#include "ROL_Types.hpp"

/** @ingroup step_group
    \class ROL::Newton_U
    \brief Provides the interface to compute optimization steps
           with Newton's method globalized using line search.
*/

namespace ROL {

template<typename Real>
class Newton_U : public DescentDirection_U<Real> {
public:
  Newton_U() {}

  void compute( Vector<Real> &s, Real &snorm, Real &sdotg, int &iter, int &flag,
          const Vector<Real> &x, const Vector<Real> &g, Objective<Real> &obj) override {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    // Compute unconstrained step
    obj.invHessVec(s,g,x,tol);
    //sdotg = -s.dot(g.dual());
    sdotg = -s.apply(g);
    if (sdotg >= static_cast<Real>(0)) {
      s.set(g.dual());
      //sdotg = -s.dot(g.dual());
      sdotg = -s.apply(g);
    }
    s.scale(static_cast<Real>(-1));
    snorm = s.norm();
    iter  = 0;
    flag  = 0;
  }

  std::string printName(void) const override {
    std::string name = "Newton's Method";
    return name;
  }
}; // class ROL::TypeU::LineSearch::Newton

} // namespace ROL

#endif
