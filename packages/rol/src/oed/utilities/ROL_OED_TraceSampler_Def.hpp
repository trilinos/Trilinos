// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_TRACESAMPLER_DEF_HPP
#define ROL_OED_TRACESAMPLER_DEF_HPP

namespace ROL {
namespace OED {

template<typename Real>
void TraceSampler<Real>::setInStorage(const Vector<Real> &g, const std::vector<Real> &param) {
  g_storage_->set(g,param);
}

template<typename Real>
bool TraceSampler<Real>::getFromStorage(Vector<Real> &g, const std::vector<Real> &param) const {
  return g_storage_->get(g,param);
}

template<typename Real>
TraceSampler<Real>::TraceSampler()
  : ProfiledClass<Real,std::string>("OED::TraceSampler") {
  g_storage_ = makePtr<SampledVector<Real>>();
}

template<typename Real>
TraceSampler<Real>::TraceSampler(const Ptr<Vector<Real>> &theta)
  : ProfiledClass<Real,std::string>("OED::TraceSampler") {
  startTimer("TraceSampler");
  Ptr<Vector<Real>> g = theta->dual().clone();
  std::vector<Real> param(1);
  g_storage_ = makePtr<SampledVector<Real>>();
  const int size = theta->dimension();
  for (int i = 0; i < size; ++i) {
    param[0] = static_cast<Real>(i);
    g_storage_->set(*g->basis(i),param);
  }
  stopTimer("TraceSampler");
}

template<typename Real>
void TraceSampler<Real>::generate(Vector<Real> &g) const {
  throw Exception::NotImplemented(">>> ROL::OED::TraceSampler::generate : Cannot generate samples from the default TraceSampler!");
}

template<typename Real>
void TraceSampler<Real>::get(Vector<Real> &F, const std::vector<Real> &param) {
  startTimer("get");
  bool isComputed = g_storage_->get(F,param);
  if (!isComputed) {
    throw Exception::NotImplemented(">>> ROL::OED::TraceSampler::get : Cannot generate samples from the default TraceSampler!");
  }
  stopTimer("get");
}

} // End OED Namespace
} // End ROL Namespace

#endif
