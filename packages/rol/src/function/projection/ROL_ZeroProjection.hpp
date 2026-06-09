// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ZEROPROJECTION_H
#define ROL_ZEROPROJECTION_H

#include "ROL_PolyhedralProjection.hpp"
#include <iostream>

namespace ROL {

template<typename Real>
class ZeroProjection : public PolyhedralProjection<Real> {
public:

  ZeroProjection() : PolyhedralProjection<Real>(nullPtr) {}

  virtual void project(Vector<Real> &x, std::ostream &stream = std::cout) override {
    x.zero();
  }

  virtual void applyJacobian(Vector<Real> &v, const Vector<Real> &x) override {
    v.zero();
  }

}; // class PolyhedralProjection_Zero
} // namespace ROL

#endif
