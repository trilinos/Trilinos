// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef POISSON_NOISE_HPP
#define POISSON_NOISE_HPP

#include "ROL_Types.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_OED_Noise.hpp"
#include <string>

template <class Real>
class Poisson_Noise : public ROL::OED::Noise<Real> {
public:
  Poisson_Noise() {}

  Real evaluate(const std::vector<Real> &param) const {
    const Real zero(0), half(0.5), one(1), five(5), s0(1e-4);
    Real norm(0);
    for (const auto x : param) norm = std::max(zero,std::abs(x-half));
    return std::exp(five*norm) - one + s0;
  }
}; // class Poisson_Noise

#endif
