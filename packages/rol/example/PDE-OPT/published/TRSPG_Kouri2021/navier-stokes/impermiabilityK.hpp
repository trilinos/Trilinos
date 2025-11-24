// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IMPERMIABILITYK_HPP
#define IMPERMIABILITYK_HPP

#include <vector>
#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"
#include "../../../TOOLS/feK.hpp"

template<class Real, class DeviceType>
class Impermiability {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;

private:
  Real alphaLo_, alphaHi_, q_;

  Real value(Real z) const {
    const Real one(1);
    return alphaHi_ + (alphaLo_-alphaHi_)*z*(one+q_)/(q_+z);
  }

  Real deriv1(Real z) const {
    const Real one(1), two(2);
    return (alphaLo_-alphaHi_)*q_*(q_+one)/std::pow(q_+z,two);
  }

  Real deriv2(Real z) const {
    const Real one(1), two(2), three(3);
    return -(alphaLo_-alphaHi_)*two*q_*(q_+one)/std::pow(q_+z,three);
  }

public:
  Impermiability(ROL::ParameterList &list) {
    alphaLo_ = list.sublist("Problem").get("Minimum Impermiability", 0.0);
    alphaHi_ = list.sublist("Problem").get("Maximum Impermiability", 1e3);
    q_       = list.sublist("Problem").get("RAMP Parameter", 1e5);
  }

  void compute(scalar_view &alpha, const scalar_view z, const scalar_view pts, int deriv=0) const {
    const int c = pts.extent_int(0);
    const int p = pts.extent_int(1);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        // Compute spatially dependent viscosity
        if (deriv==1)      alpha(i,j) = deriv1(z(i,j));
        else if (deriv==2) alpha(i,j) = deriv2(z(i,j));
        else               alpha(i,j) = value(z(i,j));
      }
    }
  }
};

#endif
