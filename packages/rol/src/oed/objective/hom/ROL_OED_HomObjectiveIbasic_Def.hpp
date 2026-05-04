// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVEIBASIC_DEF_HPP
#define ROL_OED_HOMOBJECTIVEIBASIC_DEF_HPP

namespace ROL::OED::Hom {

template<typename Real>
ObjectiveIbasic<Real>::ObjectiveIbasic( const Ptr<MomentOperator<Real>>& M,
                                        const Ptr<Factors<Real>>& F,
                                        const Ptr<SampleGenerator<Real>>& sampler,
                                        bool storage)
  : obj_(makePtr<ObjectivePV<Real>>(M,F,storage)),
    sampler_(sampler) {}

template<typename Real>
void ObjectiveIbasic<Real>::update(const Vector<Real>& z, UpdateType type, int iter) {
  obj_->update(z,type,iter);
}

template<typename Real>
Real ObjectiveIbasic<Real>::value( const Vector<Real> &z, Real &tol ) {
  Real myval(0), val(0);
  for (int i=0; i<sampler_->numMySamples(); ++i) {
    obj_->setParameter(sampler_->getMyPoint(i));
    myval += sampler_->getMyWeight(i) * obj_->value(z,tol);
  }
  sampler_->sumAll(&myval,&val,1);
  return val;
}

template<typename Real>
void ObjectiveIbasic<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (g_==nullPtr) g_ = g.clone();
  g_->zero();
  for (int i=0; i<sampler_->numMySamples(); ++i) {
    obj_->setParameter(sampler_->getMyPoint(i));
    obj_->gradient(g,z,tol);
    g_->axpy(sampler_->getMyWeight(i),g);
  }
  g.zero();
  sampler_->sumAll(*g_,g);
}

template<typename Real>
void ObjectiveIbasic<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (g_==nullPtr) g_ = hv.clone();
  g_->zero();
  for (int i=0; i<sampler_->numMySamples(); ++i) {
    obj_->setParameter(sampler_->getMyPoint(i));
    obj_->hessVec(hv,v,z,tol);
    g_->axpy(sampler_->getMyWeight(i),hv);
  }
  hv.zero();
  sampler_->sumAll(*g_,hv);
}

} // END ROL::OED::Hom Namespace

#endif
