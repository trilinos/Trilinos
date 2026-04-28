// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HETOBJECTIVEI_DEF_HPP
#define ROL_OED_HETOBJECTIVEI_DEF_HPP

#include "ROL_OED_HetObjectiveIbasic.hpp"
#include "ROL_OED_HetObjectiveItrans.hpp"

namespace ROL::OED::Het {

template<typename Real>
ObjectiveI<Real>::ObjectiveI( const Ptr<MomentOperator<Real>>& M0,
                              const Ptr<MomentOperator<Real>>& M1,
                              const Ptr<Factors<Real>>& F,
                              const Ptr<SampleGenerator<Real>>& sampler,
                              const Ptr<TraceSampler<Real>>& ts,
                              const std::vector<Real>& wt,
                              bool storage) {
  obj_ = makePtr<ObjectiveItrans<Real>>(M0,M1,F,sampler,ts,wt,storage);
}

template<typename Real>
ObjectiveI<Real>::ObjectiveI( const Ptr<MomentOperator<Real>>& M0,
                              const Ptr<MomentOperator<Real>>& M1,
                              const Ptr<Factors<Real>>& F,
                              const Ptr<SampleGenerator<Real>>& sampler,
                              bool storage, bool useBasic) {
  const unsigned nobs = F->numObservations();
  const unsigned nsamp = sampler->numMySamples();
  const unsigned dim = F->getTheta()->dimension();
  if (nobs*nsamp < dim || useBasic)
    obj_ = makePtr<ObjectiveIbasic<Real>>(M0,M1,F,sampler,storage);
  else {
    auto ts = makePtr<TraceSampler<Real>>(F->getTheta());
    std::vector<Real> wt(dim,1);
    obj_ = makePtr<ObjectiveItrans<Real>>(M0,M1,F,sampler,ts,wt,storage);
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

} // END ROL::OED::Het Namespace

#endif
