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

#ifndef PDE_GINZBURGLANDAU_EX01_HPP
#define PDE_GINZBURGLANDAU_EX01_HPP

#include "pde_ginzburg-landau.hpp"


template <class Real>
class PDE_GinzburgLandau_ex01 : public PDE_GinzburgLandau<Real> {
public:
  PDE_GinzburgLandau_ex01(Teuchos::ParameterList &parlist) : PDE_GinzburgLandau<Real>(parlist) {}

  void evaluateMagneticPotential(std::vector<Real> &Ax, const std::vector<Real> &x) const override {
    Ax[0] = static_cast<Real>(0);
    Ax[1] = static_cast<Real>(-0.5)*x[0];
  }

  Real evaluateNeumann(const std::vector<Real> &x, const int component) const override {
    const Real pi(M_PI);
    return (component==0) ? std::sin(pi*x[0])*std::cos(pi*x[1])
                          : std::cos(pi*x[0])*std::sin(pi*x[1]);
  }

  Real evaluateForce(const std::vector<Real> &x, const int component) const override {
    return static_cast<Real>(0);
  }

}; // PDE_GinzburgLandau_ex01

#endif
