// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HELMHOLTZ_NOISE_HPP
#define HELMHOLTZ_NOISE_HPP

#include "ROL_Types.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_OED_Noise.hpp"
#include <string>

template <class Real>
class Helmholtz_Noise : public ROL::OED::Noise<Real> {
public:
  Helmholtz_Noise() {}

  Real evaluate(const std::vector<Real> &param) const {
    const Real one(1), two(2), four(4), five(5);
    const Real x0(1.5), alpha(0.292893);
    const int d = param.size();
    Real val(0), linf(0), l1(0), absx(0), arg(0);
    for (int j = 0; j < d; ++j) {
      absx = std::abs(param[j]);
      linf = (linf > absx ? linf : absx);
      l1  += absx;
    }
    arg = std::max(l1/(two*(one-alpha)),linf);
    val = one + four*std::exp(five*(arg - x0));
//    Real val(s0), dot(0), loc(0), theta(0);
//    for (int i = 0; i < 8; ++i) {
//      dot = zero;
//      for (int j = 0; j < d; ++j) {
//        dot += std::pow(param[i],2);
//      }
//      val += ds*std::min(one,std::exp(twenty*(std::sqrt(dot)-x0)));
//    }
    return val;
  }
}; // class Helmholtz_Noise

#endif
