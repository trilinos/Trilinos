// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_GRADIENT_U_H
#define ROL_GRADIENT_U_H

#include "ROL_DescentDirection_U.hpp"
#include "ROL_Types.hpp"
#include "ROL_Step.hpp"

/** @ingroup step_group
    \class ROL::Gradient_U
    \brief Provides the interface to compute optimization steps
           with the gradient descent method globalized using line search.
*/

namespace ROL {

template<typename Real>
class Gradient_U : public DescentDirection_U<Real> {
public:
  Gradient_U() {}

  void compute( Vector<Real> &s, Real &snorm, Real &sdotg, int &iter, int &flag,
          const Vector<Real> &x, const Vector<Real> &g, Objective<Real> &obj) override {
    s.set(g.dual());
    s.scale(static_cast<Real>(-1));
    snorm = s.norm();
    //sdotg = s.dot(g.dual());
    sdotg = s.apply(g);
    iter  = 0;
    flag  = 0;
  }

  std::string printName(void) const override {
    std::string name = "Gradient Descent";
    return name;
  }
}; // class ROL::Gradient_U

} // namespace ROL
#endif
