// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MEANVALUEOBJECTIVE_DEF_HPP
#define ROL_MEANVALUEOBJECTIVE_DEF_HPP

namespace ROL {

template<typename Real>
MeanValueObjective<Real>::MeanValueObjective( const Ptr<Objective<Real>> &obj,
                                              const Ptr<SampleGenerator<Real>> &sampler) : obj_(obj) {
  std::vector<Real> param = computeSampleMean(sampler);
  obj_->setParameter(param);
}

template<typename Real>
void MeanValueObjective<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  obj_->update(x,type,iter);
}

template<typename Real>
void MeanValueObjective<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  obj_->update(x,flag,iter);
}

template<typename Real>
Real MeanValueObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  return obj_->value(x,tol);
}

template<typename Real>
void MeanValueObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  obj_->gradient(g,x,tol);
}

template<typename Real>
void MeanValueObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v,
        const Vector<Real> &x, Real &tol ) {
  obj_->hessVec(hv,v,x,tol);
}

template<typename Real>
void MeanValueObjective<Real>::precond( Vector<Real> &Pv, const Vector<Real> &v,
                      const Vector<Real> &x, Real &tol ) {
  obj_->precond(Pv,v,x,tol);
}

template<typename Real>
std::vector<Real> MeanValueObjective<Real>::computeSampleMean(const Ptr<SampleGenerator<Real>> &sampler) const {
  // Compute mean value of inputs and set parameter in objective
  int dim = sampler->getMyPoint(0).size(), nsamp = sampler->numMySamples();
  std::vector<Real> loc(dim), mean(dim), pt(dim);
  Real wt(0);
  for (int i = 0; i < nsamp; i++) {
    pt = sampler->getMyPoint(i);
    wt = sampler->getMyWeight(i);
    for (int j = 0; j < dim; j++) {
      loc[j] += wt*pt[j];
    }
  }
  sampler->sumAll(&loc[0],&mean[0],dim);
  return mean;
}

}

#endif
