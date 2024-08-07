// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_RADAMACHER_HPP
#define ROL_OED_RADAMACHER_HPP

#include "ROL_OED_TraceSampler.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class Radamacher : public TraceSampler<Real> {
private:
  using TraceSampler<Real>::startTimer;
  using TraceSampler<Real>::stopTimer;

public:
  Radamacher(const Ptr<Vector<Real>> &theta,int size);

  void generate(Vector<Real> &g) const;
  void get(Vector<Real> &F, const std::vector<Real> &param);

}; // class Radamacher

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_Radamacher_Def.hpp"

#endif
