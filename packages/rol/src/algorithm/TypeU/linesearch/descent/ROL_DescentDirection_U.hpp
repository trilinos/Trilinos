// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DESCENTDIRECTION_U_H
#define ROL_DESCENTDIRECTION_U_H

#include "ROL_Objective.hpp"

/** @ingroup step_group
    \class ROL::DescentDirection_U
    \brief Provides the interface to compute unconstrained optimization steps
           for line search.
*/

namespace ROL {

template <typename Real>
class DescentDirection_U {
public:
  virtual ~DescentDirection_U() {}

  virtual void initialize(const Vector<Real> &x, const Vector<Real> &g) {}

  virtual void compute( Vector<Real> &s, Real &snorm, Real &sdotg, int &iter, int &flag,
                  const Vector<Real> &x, const Vector<Real> &g, Objective<Real> &obj) = 0;

  virtual void update(const Vector<Real> &x, const Vector<Real> &s,
                      const Vector<Real> &gold, const Vector<Real> &gnew,
                      const Real snorm, const int iter) {}

  virtual std::string printName(void) const {
    std::string name = "Undefined";
    return name;
  }
}; // class DescentDirection_U
} // namespace ROL

#endif
