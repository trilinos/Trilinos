// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PROXOBJECTIVE_HPP
#define ROL_PROXOBJECTIVE_HPP

#include "ROL_Objective.hpp"

namespace ROL {

template<typename Real>
class ProxObjective : public Objective<Real> {
public:
  virtual void prox(Vector<Real> &x, Real gamma) = 0;
};

}

#endif
