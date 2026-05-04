// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVEI_DEF_HPP
#define ROL_OED_HOMOBJECTIVEI_DEF_HPP

#include "ROL_OED_HomObjectiveIbasic.hpp"
#include "ROL_OED_HomObjectiveItrans.hpp"

namespace ROL::OED::Hom {

template<typename Real>
ObjectiveI<Real>::ObjectiveI( const Ptr<MomentOperator<Real>>& M,
                              const Ptr<Factors<Real>>& F,
                              const Ptr<SampleGenerator<Real>>& sampler,
                              const Ptr<TraceSampler<Real>>& ts,
                              const std::vector<Real>& wt,
                              bool storage) {
  obj_ = makePtr<ObjectiveItrans<Real>>(M,F,sampler,ts,wt,storage);
}

template<typename Real>
ObjectiveI<Real>::ObjectiveI( const Ptr<MomentOperator<Real>>& M,
                              const Ptr<Factors<Real>>& F,
                              const Ptr<SampleGenerator<Real>>& sampler,
                              bool storage, bool useBasic) {
  const unsigned nobs = F->numObservations();
  const unsigned nsamp = sampler->numMySamples();
  const unsigned dim = F->getTheta()->dimension();
  if (nobs*nsamp < dim || useBasic)
    obj_ = makePtr<ObjectiveIbasic<Real>>(M,F,sampler,storage);
  else {
    auto ts = makePtr<TraceSampler<Real>>(F->getTheta());
    std::vector<Real> wt(dim,1);
    obj_ = makePtr<ObjectiveItrans<Real>>(M,F,sampler,ts,wt,storage);
  }
}


template<typename Real>
void ObjectiveI<Real>::update(const Vector<Real>& z, UpdateType type, int iter) {
  obj_->update(z,type,iter);
}

template<typename Real>
Real ObjectiveI<Real>::value( const Vector<Real> &z, Real &tol ) {
  return obj_->value(z,tol);
}

template<typename Real>
void ObjectiveI<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  obj_->gradient(g,z,tol);
}

template<typename Real>
void ObjectiveI<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  obj_->hessVec(hv,v,z,tol);
}

} // END ROL::OED::Hom Namespace

#endif
