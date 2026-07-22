// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ZEROPROXOBJECTIVE_HPP
#define ROL_ZEROPROXOBJECTIVE_HPP

#include "ROL_ProxObjective.hpp"

namespace ROL {

template<typename Real>
class ZeroProxObjective : public ProxObjective<Real> {
public:

  Real value(const Vector<Real> x, Real &tol) override {
    return static_cast<Real>(0);
  }

  void gradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) override {
    g.zero();
  }

  void prox(Vector<Real> &x, Real gamma) override {}

};

}

#endif
