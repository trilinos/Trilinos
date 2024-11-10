// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_RADAMACHER_DEF_HPP
#define ROL_OED_RADAMACHER_DEF_HPP

namespace ROL {
namespace OED {

template<typename Real>
Radamacher<Real>::Radamacher(const Ptr<Vector<Real>> &theta, int size) {
  ProfiledClass<Real,std::string>::rename("OED::Radamacher");
  startTimer("Radamacher");
  Ptr<Vector<Real>> g = theta->dual().clone();
  std::vector<Real> param(1);
  for (int i = 0; i < size; ++i) {
    param[0] = static_cast<Real>(i);
    generate(*g);
    TraceSampler<Real>::setInStorage(*g,param);
  }
  stopTimer("Radamacher");
}

template<typename Real>
void Radamacher<Real>::generate(Vector<Real> &g) const {
  startTimer("generate");
  g.randomize();
  g.applyUnary(Elementwise::Round<Real>());
  g.scale(static_cast<Real>(2));
  g.applyUnary(Elementwise::Shift<Real>(static_cast<Real>(-1)));
  stopTimer("generate");
}

template<typename Real>
void Radamacher<Real>::get(Vector<Real> &F, const std::vector<Real> &param) {
  startTimer("get");
  bool isComputed = TraceSampler<Real>::getFromStorage(F,param);
  if (!isComputed) {
    generate(F);
    TraceSampler<Real>::setInStorage(F,param);
  }
  stopTimer("get");
}

} // End OED Namespace
} // End ROL Namespace

#endif
