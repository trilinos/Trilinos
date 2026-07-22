// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pde_ginzburg-landau.hpp
    \brief Implements the local PDE interface for the optimal control of
           simplified Ginzburg-Landau.
*/

#ifndef PDE_GINZBURGLANDAU_EX02_HPP
#define PDE_GINZBURGLANDAU_EX02_HPP

#include "pde_ginzburg-landau.hpp"


template <class Real>
class PDE_GinzburgLandau_ex02 : public PDE_GinzburgLandau<Real> {
public:
  PDE_GinzburgLandau_ex02(Teuchos::ParameterList &parlist) : PDE_GinzburgLandau<Real>(parlist) {}

  void evaluateMagneticPotential(std::vector<Real> &Ax, const std::vector<Real> &x) const override {
    const Real pi(M_PI);
    Ax[0] =  std::sin(pi*x[0])*std::cos(pi*x[1]);
    Ax[1] = -std::cos(pi*x[0])*std::sin(pi*x[1]);
  }

  Real evaluateNeumann(const std::vector<Real> &x, const int component) const override {
    return static_cast<Real>(0);
  }

  Real evaluateForce(const std::vector<Real> &x, const int component) const override {
    const Real pi(M_PI), one(1), two(2);
    const Real cx = std::cos(pi*x[0]), sx = std::sin(pi*x[0]);
    const Real cy = std::cos(pi*x[1]), sy = std::sin(pi*x[1]);
    const Real pi2 = pi*pi, cx2 = cx*cx, cy2 = cy*cy, sx2 = sx*sx, sy2 = sy*sy;
    const Real rhs = pi2 + sx2*cy2 + cx2*sy2 - one + cx2 + cy2;
    return (component == 0) ? cx*(rhs - two*pi*sy2) : cy*(rhs - two*pi*sx2);
  }

}; // PDE_GinzburgLandau_ex01

#endif
