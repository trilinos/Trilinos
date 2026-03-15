// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_ROBUSTOBJECTIVE_HPP
#define ROL_OED_ROBUSTOBJECTIVE_HPP

#include "ROL_Objective.hpp"
#include "ROL_OED_DesignVector.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class RobustObjective : public Objective<Real> {
private:

public:
  RobustObjective() {}

  Real value( const Vector<Real> &x, Real &tol ) override {
    const DesignVector<Real> &xs = static_cast<const DesignVector<Real>&>(x);
    return xs.getValue();
  }
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override {
    g.zero();
    DesignVector<Real> &gs = static_cast<DesignVector<Real>&>(g);
    gs.setValue(static_cast<Real>(1));
  }
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override {
    hv.zero();
  }

}; // class RobustObjective

} // End OED Namespace
} // End ROL Namespace

#endif
