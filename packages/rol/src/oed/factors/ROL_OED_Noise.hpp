// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_NOISE_HPP
#define ROL_OED_NOISE_HPP

#include <vector>

namespace ROL {
namespace OED {

template<typename Real>
class Noise {
public:
  virtual ~Noise() {}

  // Standard deviation of noise
  virtual Real evaluate(const std::vector<Real> &param) const {
    return static_cast<Real>(1);
  }

  // Moment generating function of noise evaluated at t=2
  virtual Real mgf2(const std::vector<Real> &param) const {
    return static_cast<Real>(2);
  }
}; // class Noise

} // End OED Namespace
} // End ROL Namespace

#endif
