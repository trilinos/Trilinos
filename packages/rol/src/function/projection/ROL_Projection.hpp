// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PROJECTION_H
#define ROL_PROJECTION_H

#include <iostream>
#include "ROL_Vector.hpp"

namespace ROL {

template<typename Real>
class Projection {
public:

  virtual ~Projection() {}

  virtual void project(Vector<Real> &x, std::ostream &stream = std::cout) = 0;

  virtual void applyJacobian(Vector<Real> &v, const Vector<Real> &x) {
    throw Exception::NotImplemented(">>> ROL::Projection::applyJacobian: Not Implemented!");
  }

}; // class Projection

} // namespace ROL

#endif
